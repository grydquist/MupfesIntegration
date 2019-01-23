!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This is for solving Saint Venant-Kirchhoff (svk) equations. 
!     For usage and testing see TEST_SVKMOD at the end of this file.
!      
!---------------------------------------------------------------------
      MODULE SVKMOD
      USE EQMOD
      IMPLICIT NONE

!     Is used to initialize equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 1,
     2   itgO = 2,
     3   lsAl = LS_TYPE_CG,
     4   sym  = "SV")

!     Equation type. Extending solution control 
      TYPE, EXTENDS(eqType) :: svkType
!        Displacement, velocity, ...
         TYPE(gVarType), POINTER :: U => NULL()
      CONTAINS
!        Setups all structure
         PROCEDURE :: setup => setupSvk
!        Evaulate the governing equation
         PROCEDURE :: eval3 => eval3Svk
         PROCEDURE :: eval2 => eval2Svk
!        Overridden procedures (no implementation)
         PROCEDURE :: bEval => bEvalSvk
      END TYPE svkType

      INTERFACE svkType
         PROCEDURE :: newSvk, newSvkFl
      END INTERFACE svkType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newSvk(dmn, sct, mName, nBc, ls, itg) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: mName
      INTEGER, INTENT(IN) :: nBc
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      TYPE(itgType), INTENT(IN), OPTIONAL :: itg
      TYPE(svkType) :: eq
      
      CALL eq%new(eqSp, dmn, sct, mName, nBc, ls, itg)

      RETURN
      END FUNCTION newSvk
!---------------------------------------------------------------------
      FUNCTION newSvkFl(dmn, lst) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(svkType) :: eq
      
      CALL eq%new(eqSp, dmn, lst)

      RETURN
      END FUNCTION newSvkFl
!---------------------------------------------------------------------
      SUBROUTINE setupSvk(eq, var)
      CLASS(svkType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)

      IF (PRESENT(var)) THEN
         IF (SIZE(var) .NE. 1) io%e = "setupLed: Invalid var size"
         IF (var(1)%dof .NE. nsd) io%e = "setupLed: Invalid var nU"
         eq%var => var
      ELSE
         eq%var(1) = gVarType(nsd,'Velocity',eq%dmn)
      END IF
      eq%U => eq%var(1)
      eq%conf = eqConf_ref

      IF (eq%mat%rho() .LE. 0D0) io%e = "setupSvk: rho <= 0"
      IF (eq%mat%E() .LE. 0D0)   io%e = "setupSvk: E <= 0"
      IF (eq%mat%nu() .GE. 5D-1) io%e = "setupSvk: nu >= 0.5"
      IF (eq%mat%nu() .LE. 0D0)  io%e = "setupSvk: nu <= 0.0"

      RETURN
      END SUBROUTINE setupSvk
!---------------------------------------------------------------------
      PURE SUBROUTINE eval3Svk(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(svkType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
      
      INTEGER a, b, Ac
      REAL(KIND=8) NxdNx, lambda, mu, rho, elM, nu, C(nsd,nsd), lDm,
     2   S(nsd,nsd), E(nsd,nsd), T1, kappa, trE, fs(nsd,nsd), 
     3   NxF(nsd,eNoN), NxSNx, FFt(nsd,nsd), amd, wl, F(nsd,nsd), 
     4   ud(nsd), bf(nsd), dmp

      rho   = eq%mat%rho()
      elM   = eq%mat%E()
      nu    = eq%mat%nu() 
      bf(1) = eq%mat%fx()
      bf(2) = eq%mat%fy()
      bf(3) = eq%mat%fz()
      dmp   = eq%mat%eta()/rho
      
      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      kappa  = lambda + 2D0/3D0*mu
      lDm    = lambda/mu
      T1     = eq%itg%af*eq%itg%beta*dt*dt
      amd    = (eq%itg%am + eq%itg%af*eq%itg%gam*dt*dmp)/T1*rho
      wl     = w*T1*mu

      f      = 0D0
      f(1,1) = 1D0
      f(2,2) = 1D0
      f(3,3) = 1D0
      ud     = -bf
      DO a=1, eNoN
         Ac = eqN(a)
         ud(1)  = ud(1) + N(a)*(eq%U%A%v(1,Ac) + dmp*eq%U%v(1,Ac))
         ud(2)  = ud(2) + N(a)*(eq%U%A%v(2,Ac) + dmp*eq%U%v(2,Ac))
         ud(3)  = ud(3) + N(a)*(eq%U%A%v(3,Ac) + dmp*eq%U%v(3,Ac))

         f(1,1) = f(1,1) + Nx(1,a)*eq%U%D%v(1,Ac)
         f(2,1) = f(2,1) + Nx(1,a)*eq%U%D%v(2,Ac)
         f(3,1) = f(3,1) + Nx(1,a)*eq%U%D%v(3,Ac)
         f(1,2) = f(1,2) + Nx(2,a)*eq%U%D%v(1,Ac)
         f(2,2) = f(2,2) + Nx(2,a)*eq%U%D%v(2,Ac)
         f(3,2) = f(3,2) + Nx(2,a)*eq%U%D%v(3,Ac)
         f(1,3) = f(1,3) + Nx(3,a)*eq%U%D%v(1,Ac)
         f(2,3) = f(2,3) + Nx(3,a)*eq%U%D%v(2,Ac)
         f(3,3) = f(3,3) + Nx(3,a)*eq%U%D%v(3,Ac)
      END DO

      C(1,1) = F(1,1)*F(1,1) + F(2,1)*F(2,1) + F(3,1)*F(3,1)
      C(1,2) = F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,1)*F(3,2)
      C(1,3) = F(1,1)*F(1,3) + F(2,1)*F(2,3) + F(3,1)*F(3,3)
      C(2,1) = C(1,2)
      C(2,2) = F(1,2)*F(1,2) + F(2,2)*F(2,2) + F(3,2)*F(3,2)
      C(2,3) = F(1,2)*F(1,3) + F(2,2)*F(2,3) + F(3,2)*F(3,3)
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)
      C(3,3) = F(1,3)*F(1,3) + F(2,3)*F(2,3) + F(3,3)*F(3,3)

      E(1,1) = 5D-1*C(1,1) - 5D-1
      E(1,2) = 5D-1*C(1,2)
      E(1,3) = 5D-1*C(1,3)
      E(2,1) = 5D-1*C(2,1)
      E(2,2) = 5D-1*C(2,2) - 5D-1
      E(2,3) = 5D-1*C(2,3)
      E(3,1) = 5D-1*C(3,1)
      E(3,2) = 5D-1*C(3,2)
      E(3,3) = 5D-1*C(3,3) - 5D-1
     
      trE = E(1,1) + E(2,2) + E(3,3)

      S(1,1) = lambda*trE + 2D0*mu*E(1,1)
      S(1,2) = 2D0*mu*E(1,2)
      S(1,3) = 2D0*mu*E(1,3)
      S(2,1) = S(1,2)
      S(2,2) = lambda*trE + 2D0*mu*E(2,2)
      S(2,3) = 2D0*mu*E(2,3)
      S(3,1) = S(1,3)
      S(3,2) = S(2,3)
      S(3,3) = lambda*trE + 2D0*mu*E(3,3)

      FFt(1,1) = F(1,1)*F(1,1) + F(1,2)*F(1,2) + F(1,3)*F(1,3)
      FFt(1,2) = F(1,1)*F(2,1) + F(1,2)*F(2,2) + F(1,3)*F(2,3)
      FFt(1,3) = F(1,1)*F(3,1) + F(1,2)*F(3,2) + F(1,3)*F(3,3)
      FFt(2,1) = FFt(1,2)
      FFt(2,2) = F(2,1)*F(2,1) + F(2,2)*F(2,2) + F(2,3)*F(2,3)
      FFt(2,3) = F(2,1)*F(3,1) + F(2,2)*F(3,2) + F(2,3)*F(3,3)
      FFt(3,1) = FFt(1,3)
      FFt(3,2) = FFt(2,3)
      FFt(3,3) = F(3,1)*F(3,1) + F(3,2)*F(3,2) + F(3,3)*F(3,3)
      
      fs(1,1) = F(1,1)*S(1,1) + F(1,2)*S(2,1) + F(1,3)*S(3,1)
      fs(1,2) = F(1,1)*S(1,2) + F(1,2)*S(2,2) + F(1,3)*S(3,2)
      fs(1,3) = F(1,1)*S(1,3) + F(1,2)*S(2,3) + F(1,3)*S(3,3)
      fs(2,1) = F(2,1)*S(1,1) + F(2,2)*S(2,1) + F(2,3)*S(3,1)
      fs(2,2) = F(2,1)*S(1,2) + F(2,2)*S(2,2) + F(2,3)*S(3,2)
      fs(2,3) = F(2,1)*S(1,3) + F(2,2)*S(2,3) + F(2,3)*S(3,3)
      fs(3,1) = F(3,1)*S(1,1) + F(3,2)*S(2,1) + F(3,3)*S(3,1)
      fs(3,2) = F(3,1)*S(1,2) + F(3,2)*S(2,2) + F(3,3)*S(3,2)
      fs(3,3) = F(3,1)*S(1,3) + F(3,2)*S(2,3) + F(3,3)*S(3,3)
      
      DO a=1, eNoN
         NxF(1,a) = F(1,1)*Nx(1,a) + F(1,2)*Nx(2,a) + F(1,3)*Nx(3,a)
         NxF(2,a) = F(2,1)*Nx(1,a) + F(2,2)*Nx(2,a) + F(2,3)*Nx(3,a)
         NxF(3,a) = F(3,1)*Nx(1,a) + F(3,2)*Nx(2,a) + F(3,3)*Nx(3,a)
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1)
     2      + Nx(1,a)*fs(1,1) + Nx(2,a)*fs(1,2) + Nx(3,a)*fs(1,3))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2)
     2      + Nx(1,a)*fs(2,1) + Nx(2,a)*fs(2,2) + Nx(3,a)*fs(2,3))
         
         lR(3,a) = lR(3,a) + w*(rho*N(a)*ud(3)
     2      + Nx(1,a)*fs(3,1) + Nx(2,a)*fs(3,2) + Nx(3,a)*fs(3,3))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)
            NxSNx = Nx(1,a)*S(1,1)*Nx(1,b) + Nx(1,a)*S(1,2)*Nx(2,b) 
     2            + Nx(1,a)*S(1,3)*Nx(3,b) + Nx(2,a)*S(2,1)*Nx(1,b) 
     3            + Nx(2,a)*S(2,2)*Nx(2,b) + Nx(2,a)*S(2,3)*Nx(3,b) 
     4            + Nx(3,a)*S(3,1)*Nx(1,b) + Nx(3,a)*S(3,2)*Nx(2,b) 
     5            + Nx(3,a)*S(3,3)*Nx(3,b)
            
            T1    = (amd*N(a)*N(b) + NxSNx)/mu
            
            lK(1,a,b) = lK(1,a,b) + wl*(T1 + FFt(1,1)*NxdNx
     2         + (1D0 + lDm)*NxF(1,b)*NxF(1,a))

            lK(2,a,b) = lK(2,a,b) + wl*(NxF(1,b)*NxF(2,a) 
     2         + FFt(1,2)*NxdNx + lDm*NxF(1,a)*NxF(2,b))

            lK(3,a,b) = lK(3,a,b) + wl*(NxF(1,b)*NxF(3,a) 
     2         + FFt(1,3)*NxdNx + lDm*NxF(1,a)*NxF(3,b))

            lK(4,a,b) = lK(4,a,b) + wl*(NxF(2,b)*NxF(1,a) 
     2         + FFt(2,1)*NxdNx + lDm*NxF(2,a)*NxF(1,b))

            lK(5,a,b) = lK(5,a,b) + wl*(T1 + FFt(2,2)*NxdNx
     2         + (1D0 + lDm)*NxF(2,b)*NxF(2,a))
            
            lK(6,a,b) = lK(6,a,b) + wl*(NxF(2,b)*NxF(3,a) 
     2         + FFt(2,3)*NxdNx + lDm*NxF(2,a)*NxF(3,b))

            lK(7,a,b) = lK(7,a,b) + wl*(NxF(3,b)*NxF(1,a)
     2         + FFt(3,1)*NxdNx + lDm*NxF(3,a)*NxF(1,b))

            lK(8,a,b) = lK(8,a,b) + wl*(NxF(3,b)*NxF(2,a)
     2         + FFt(3,2)*NxdNx + lDm*NxF(3,a)*NxF(2,b))
            
            lK(9,a,b) = lK(9,a,b) + wl*(T1 + FFt(3,3)*NxdNx
     2         + (1D0 + lDm)*NxF(3,b)*NxF(3,a))
         END DO
      END DO

      RETURN
      END SUBROUTINE eval3Svk
!---------------------------------------------------------------------
      PURE SUBROUTINE eval2Svk(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(svkType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
      
      INTEGER a, b, Ac
      REAL(KIND=8) NxdNx, lambda, mu, rho, elM, nu, C(nsd,nsd), lDm,
     2   S(nsd,nsd), E(nsd,nsd), T1, kappa, trE, fs(nsd,nsd), 
     3   NxF(nsd,eNoN), NxSNx, FFt(nsd,nsd), amd, wl, F(nsd,nsd), 
     4   ud(nsd), bf(nsd), dmp

      rho   = eq%mat%rho()
      elM   = eq%mat%E()
      nu    = eq%mat%nu() 
      bf(1) = eq%mat%fx()
      bf(2) = eq%mat%fy()
      dmp   = eq%mat%eta()/rho
      
      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      kappa  = lambda + 2D0/3D0*mu
      lDm    = lambda/mu
      T1     = eq%itg%af*eq%itg%beta*dt*dt
      amd    = (eq%itg%am + eq%itg%af*eq%itg%gam*dt*dmp)/T1*rho
      wl     = w*T1*mu

      f      = 0D0
      f(1,1) = 1D0
      f(2,2) = 1D0
      ud     = -bf
      DO a=1, eNoN
         Ac = eqN(a)
         ud(1)  = ud(1) + N(a)*(eq%U%A%v(1,Ac) + dmp*eq%U%v(1,Ac))
         ud(2)  = ud(2) + N(a)*(eq%U%A%v(2,Ac) + dmp*eq%U%v(2,Ac))

         f(1,1) = f(1,1) + Nx(1,a)*eq%U%D%v(1,Ac)
         f(2,1) = f(2,1) + Nx(1,a)*eq%U%D%v(2,Ac)
         f(1,2) = f(1,2) + Nx(2,a)*eq%U%D%v(1,Ac)
         f(2,2) = f(2,2) + Nx(2,a)*eq%U%D%v(2,Ac)
      END DO

      C(1,1) = F(1,1)*F(1,1) + F(2,1)*F(2,1)
      C(1,2) = F(1,1)*F(1,2) + F(2,1)*F(2,2)
      C(2,1) = C(1,2)
      C(2,2) = F(1,2)*F(1,2) + F(2,2)*F(2,2)

      E(1,1) = 5D-1*C(1,1) - 5D-1
      E(1,2) = 5D-1*C(1,2)
      E(2,1) = 5D-1*C(2,1)
      E(2,2) = 5D-1*C(2,2) - 5D-1
     
      trE = E(1,1) + E(2,2)

      S(1,1) = lambda*trE + 2D0*mu*E(1,1)
      S(1,2) = 2D0*mu*E(1,2)
      S(2,1) = S(1,2)
      S(2,2) = lambda*trE + 2D0*mu*E(2,2)

      FFt(1,1) = F(1,1)*F(1,1) + F(1,2)*F(1,2)
      FFt(1,2) = F(1,1)*F(2,1) + F(1,2)*F(2,2)
      FFt(2,1) = FFt(1,2)
      FFt(2,2) = F(2,1)*F(2,1) + F(2,2)*F(2,2)
      
      fs(1,1) = F(1,1)*S(1,1) + F(1,2)*S(2,1)
      fs(1,2) = F(1,1)*S(1,2) + F(1,2)*S(2,2)
      fs(2,1) = F(2,1)*S(1,1) + F(2,2)*S(2,1)
      fs(2,2) = F(2,1)*S(1,2) + F(2,2)*S(2,2)
      
      DO a=1, eNoN
         NxF(1,a) = F(1,1)*Nx(1,a) + F(1,2)*Nx(2,a)
         NxF(2,a) = F(2,1)*Nx(1,a) + F(2,2)*Nx(2,a)
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(1,a)*fs(1,1) 
     2      + Nx(2,a)*fs(1,2))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + Nx(1,a)*fs(2,1) 
     2      + Nx(2,a)*fs(2,2))
         
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)
            NxSNx = Nx(1,a)*S(1,1)*Nx(1,b) + Nx(1,a)*S(1,2)*Nx(2,b) 
     2            + Nx(2,a)*S(2,1)*Nx(1,b) + Nx(2,a)*S(2,2)*Nx(2,b)
            
            T1    = (amd*N(a)*N(b) + NxSNx)/mu

            lK(1,a,b) = lK(1,a,b) + wl*(T1 + FFt(1,1)*NxdNx
     2         + (1D0 + lDm)*NxF(1,b)*NxF(1,a))

            lK(2,a,b) = lK(2,a,b) + wl*(NxF(1,b)*NxF(2,a) 
     2         + FFt(1,2)*NxdNx + lDm*NxF(1,a)*NxF(2,b))

            lK(3,a,b) = lK(3,a,b) + wl*(NxF(2,b)*NxF(1,a) 
     2         + FFt(2,1)*NxdNx + lDm*NxF(2,a)*NxF(1,b))

            lK(4,a,b) = lK(4,a,b) + wl*(T1 + FFt(2,2)*NxdNx
     2         + (1D0 + lDm)*NxF(2,b)*NxF(2,a))
         END DO
      END DO

      RETURN
      END SUBROUTINE eval2Svk
!---------------------------------------------------------------------
      PURE SUBROUTINE bEvalSvk(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(svkType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)


      RETURN
      END SUBROUTINE bEvalSvk

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_SVKMOD(dmn, sct, gg, res)
      TYPE(dmnType), INTENT(INOUT) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(varType), INTENT(INOUT) :: res

      INTEGER bType
      TYPE(svkType) eq

      cTS   = 0 
      eq    = svkType(dmn, sct, 'Steel', 1) ! ., ., mat_name, nBc
      bType = IBSET(0,bType_Dir)
      bType = IBSET(bType,bType_impD)
      eq%bc(1) = bcType(dmn, "Face_1", gg, bType) ! ., fa_name, ., bType
      CALL eq%ini()
      DO cTS=1, nTS
         DO WHILE(.NOT.eq%satisfied())
            CALL eq%solve()
         END DO
         CALL eq%update()
         CALL eq%save()
      END DO
      res%s = eq%U%v(2,:)
      CALL eq%free()

      RETURN
      END SUBROUTINE TEST_SVKMOD
      END MODULE SVKMOD
