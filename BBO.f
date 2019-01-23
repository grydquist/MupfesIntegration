!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This is for solving Basset–Boussinesq–Oseen (BBO) equations. 
!     For usage and testing see TEST_BBOMOD at the end of this file.
!      
!---------------------------------------------------------------------
      MODULE BBOMOD
      USE INSMOD
      IMPLICIT NONE

!     Is used to initialize equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 1,
     2   itgO = 1,
     3   lsAl = LS_TYPE_GMRES,
     5   sym  = "BB")

!     Equation type. Extending solution control 
      TYPE, EXTENDS(eqType) :: bboType
!        velocity
         TYPE(gVarType), POINTER :: U => NULL()
!        Pointer to the fluid velocity
         TYPE(gVarType), POINTER :: Uf => NULL()
!        Pointer to the fluid pressure field
         TYPE(varType), POINTER :: P => NULL()
!        Fluid properties
         TYPE(matType), POINTER, PUBLIC :: fmat => NULL()
      CONTAINS
!        Setups all structure
         PROCEDURE :: setup => setupBbo
!        Evaulate the governing equation
         PROCEDURE :: eval3 => eval3Bbo
         PROCEDURE :: eval2 => eval2Bbo
!        Overridden procedures (no implementation)
         PROCEDURE :: bEval => bEvalBbo
      END TYPE bboType

      INTERFACE bboType
         PROCEDURE :: newBbo, newBboFl
      END INTERFACE bboType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newBbo(ns, sct, mName, nBc, ls) RESULT(eq)
      TYPE(insType), INTENT(IN), TARGET :: ns
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: mName
      INTEGER, INTENT(IN) :: nBc
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      TYPE(bboType) :: eq
      
      CALL eq%new(eqSp, ns%dmn, sct, mName, nBc, ls)
      eq%Uf   => ns%U
      eq%P    => ns%P%varType
      eq%fmat => ns%mat

      RETURN
      END FUNCTION newBbo
!---------------------------------------------------------------------
      FUNCTION newBboFl(ns, lst) RESULT(eq)
      TYPE(insType), INTENT(IN), TARGET :: ns
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(bboType) :: eq
      
      CALL eq%new(eqSp, ns%dmn, lst)
      eq%Uf   => ns%U
      eq%P    => ns%P%varType
      eq%fmat => ns%mat

      RETURN
      END FUNCTION newBboFl
!---------------------------------------------------------------------
      SUBROUTINE setupBbo(eq, var)
      CLASS(bboType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
      
      IF (PRESENT(var)) THEN
         IF (SIZE(var) .NE. 1) io%e = "setupBbo: Invalid var size"
         IF (var(1)%dof .NE. nsd) io%e = "setupBbo: Invalid var nU"
         eq%var => var
      ELSE
         eq%var(1) = gVarType(nsd,'particle_velocity',eq%dmn)
      END IF
      eq%U => eq%var(1)

      IF (eq%mat%rho() .LE. 0D0) io%e = "setupBbo: rho <= 0"
      IF (eq%mat%D() .LE. 0D0)   io%e = "setupBbo: D <= 0"

      RETURN
      END SUBROUTINE setupBbo
!---------------------------------------------------------------------
      PURE SUBROUTINE eval3Bbo(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(bboType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
 
      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,1D0,4D0/)
      INTEGER a, b, Ac
      REAL(KIND=8) wr, T1, T2, amd, wl, udNx(eNoN), res(nsd), tauM, kU, 
     2   nu, tauB, up(nsd), updNx(eNoN), tmpV(nsd), NxdNx, rM(nsd,nsd), 
     3   tauMs, wlnu, ud(nsd), u(nsd), ux(nsd,nsd), bf(nsd), st, 
     4   px(nsd), rho, pRho, pDia
 
      pDia  = eq%mat%D()
      pRho  = eq%mat%rho()
      rho   = eq%fmat%rho()
      nu    = eq%fmat%mu()/rho
      bf(1) = eq%mat%fx() - eq%fmat%fx()
      bf(2) = eq%mat%fy() - eq%fmat%fy()
      bf(3) = eq%mat%fz() - eq%fmat%fz()
      st    = 18D0*nu/pDia/pDia/(pRho + 5D-1*rho)
      
      T1    = eq%itg%af*eq%itg%gam*dt
      amd   = eq%itg%am/T1
      wr    = w*(pRho + 5D-1*rho)
      wl    = wr*T1
      
      ud = 0D0
      u  = 0D0
      up = 0D0
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         up(1) = up(1) + N(a)*eq%U%v(1,Ac)
         up(2) = up(2) + N(a)*eq%U%v(2,Ac)
         up(3) = up(3) + N(a)*eq%U%v(3,Ac)

         ud(1) = ud(1) + N(a)*eq%Uf%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%Uf%A%v(2,Ac)
         ud(3) = ud(3) + N(a)*eq%Uf%A%v(3,Ac)
 
         u(1) = u(1) + N(a)*eq%Uf%v(1,Ac)
         u(2) = u(2) + N(a)*eq%Uf%v(2,Ac)
         u(3) = u(3) + N(a)*eq%Uf%v(3,Ac)
 
         px(1) = px(1) + Nx(1,a)*eq%P%s(Ac)
         px(2) = px(2) + Nx(2,a)*eq%P%s(Ac)
         px(3) = px(3) + Nx(3,a)*eq%P%s(Ac)

         ux(1,1) = ux(1,1) + Nx(1,a)*eq%Uf%v(1,Ac)
         ux(2,1) = ux(2,1) + Nx(2,a)*eq%Uf%v(1,Ac)
         ux(3,1) = ux(3,1) + Nx(3,a)*eq%Uf%v(1,Ac)
         ux(1,2) = ux(1,2) + Nx(1,a)*eq%Uf%v(2,Ac)
         ux(2,2) = ux(2,2) + Nx(2,a)*eq%Uf%v(2,Ac)
         ux(3,2) = ux(3,2) + Nx(3,a)*eq%Uf%v(2,Ac)
         ux(1,3) = ux(1,3) + Nx(1,a)*eq%Uf%v(3,Ac)
         ux(2,3) = ux(2,3) + Nx(2,a)*eq%Uf%v(3,Ac)
         ux(3,3) = ux(3,3) + Nx(3,a)*eq%Uf%v(3,Ac)
      END DO
 
      bf(1) = st*u(1) + (-rho*px(1) + 5D-1*rho*(ud(1)
     2   + up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
     3   + bf(1))/(pRho + 5D-1*rho)
      bf(2) = st*u(2) + (-rho*px(2) + 5D-1*rho*(ud(2)
     2   + up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
     3   + bf(2))/(pRho + 5D-1*rho)
      bf(3) = st*u(3) + (-rho*px(3) + 5D-1*rho*(ud(3)
     2   + up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))
     3   + bf(3))/(pRho + 5D-1*rho)
       
      u  = up
      ud = 0D0
      ux = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         ud(1) = ud(1) + N(a)*eq%U%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%U%A%v(2,Ac)
         ud(3) = ud(3) + N(a)*eq%U%A%v(3,Ac)
 
         ux(1,1) = ux(1,1) + Nx(1,a)*eq%U%v(1,Ac)
         ux(2,1) = ux(2,1) + Nx(2,a)*eq%U%v(1,Ac)
         ux(3,1) = ux(3,1) + Nx(3,a)*eq%U%v(1,Ac)
         ux(1,2) = ux(1,2) + Nx(1,a)*eq%U%v(2,Ac)
         ux(2,2) = ux(2,2) + Nx(2,a)*eq%U%v(2,Ac)
         ux(3,2) = ux(3,2) + Nx(3,a)*eq%U%v(2,Ac)
         ux(1,3) = ux(1,3) + Nx(1,a)*eq%U%v(3,Ac)
         ux(2,3) = ux(2,3) + Nx(2,a)*eq%U%v(3,Ac)
         ux(3,3) = ux(3,3) + Nx(3,a)*eq%U%v(3,Ac)
      END DO

      kU = u(1)*u(1)*ks(1,1) + u(2)*u(1)*ks(2,1) + u(3)*u(1)*ks(3,1) 
     2   + u(1)*u(2)*ks(1,2) + u(2)*u(2)*ks(2,2) + u(3)*u(2)*ks(3,2) 
     4   + u(1)*u(3)*ks(1,3) + u(2)*u(3)*ks(2,3) + u(3)*u(3)*ks(3,3)

      tauM = ct(3)/SQRT(ct(1)/dt/dt + ct(2)*kU + st*st)
      tauMs = ct(3)/SQRT(ct(2)*kU + st*st)
      
      res(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + u(3)*ux(3,1)
     2   + st*u(1) - bf(1)
      res(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + u(3)*ux(3,2) 
     2   + st*u(2) - bf(2)
      res(3) = ud(3) + u(1)*ux(1,3) + u(2)*ux(2,3) + u(3)*ux(3,3) 
     2   + st*u(3) - bf(3)

!     Discontinuity capturing nu (wx , ux)
      tmpV(1) = (u(1)*ux(1,1) + u(2)*ux(1,2) + u(3)*ux(1,3))
      tmpV(2) = (u(1)*ux(2,1) + u(2)*ux(2,2) + u(3)*ux(2,3))
      tmpV(3) = (u(1)*ux(3,1) + u(2)*ux(3,2) + u(3)*ux(3,3))
      
      nu = tmpV(1)*tmpV(1)*ks(1,1) + tmpV(2)*tmpV(1)*ks(2,1) 
     2   + tmpV(3)*tmpV(1)*ks(3,1) + tmpV(1)*tmpV(2)*ks(1,2) 
     3   + tmpV(2)*tmpV(2)*ks(2,2) + tmpV(3)*tmpV(2)*ks(3,2) 
     4   + tmpV(1)*tmpV(3)*ks(1,3) + tmpV(2)*tmpV(3)*ks(2,3) 
     5   + tmpV(3)*tmpV(3)*ks(3,3)
      IF (ISZERO(nu)) nu = eps
      nu = (u(1)*u(1) + u(2)*u(2) + u(3)*u(3))/nu
      nu = SQRT(nu*(res(1)*res(1) + res(2)*res(2) + res(3)*res(3)))/2D0
      IF (ISZERO(kU)) kU = eps
      T1 = (u(1)*u(1) + u(2)*u(2) + u(3)*u(3))/SQRT(kU)/2D0
      IF (nu .GT. T1) nu = T1

!     tauB terms tauB (up.wx , up.ux)
      up = -tauM*res
      tauB = up(1)*up(1)*ks(1,1) + up(2)*up(1)*ks(2,1) 
     2     + up(3)*up(1)*ks(3,1) + up(1)*up(2)*ks(1,2) 
     3     + up(2)*up(2)*ks(2,2) + up(3)*up(2)*ks(3,2) 
     4     + up(1)*up(3)*ks(1,3) + up(2)*up(3)*ks(2,3) 
     5     + up(3)*up(3)*ks(3,3)
      IF (ISZERO(tauB)) tauB = eps
      tauB = ct(3)/SQRT(tauB)
      
      rM(1,1) = nu*(ux(1,1) + ux(1,1)) 
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) 
      rM(3,1) = nu*(ux(3,1) + ux(1,3)) 
      rM(1,2) = nu*(ux(1,2) + ux(2,1))
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) 
      rM(3,2) = nu*(ux(3,2) + ux(2,3)) 
      rM(1,3) = nu*(ux(1,3) + ux(3,1))
      rM(2,3) = nu*(ux(2,3) + ux(3,2)) 
      rM(3,3) = nu*(ux(3,3) + ux(3,3)) 

      tmpV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
      tmpV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
      tmpV(3) = tauB*(up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))
      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)  + u(3)*Nx(3,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a) + up(3)*Nx(3,a)
         T1       = N(a) + udNx(a)*tauM

         lR(1,a) = lR(1,a) + wr*(res(1)*T1 + updNx(a)*tmpV(1) 
     2           + rM(1,1)*Nx(1,a) + rM(1,2)*Nx(2,a) + rM(1,3)*Nx(3,a))
         lR(2,a) = lR(2,a) + wr*(res(2)*T1 + updNx(a)*tmpV(2)
     2           + rM(2,1)*Nx(1,a) + rM(2,2)*Nx(2,a) + rM(2,3)*Nx(3,a))
         lR(3,a) = lR(3,a) + wr*(res(3)*T1 + updNx(a)*tmpV(3)
     2           + rM(3,1)*Nx(1,a) + rM(3,2)*Nx(2,a) + rM(3,3)*Nx(3,a))
      END DO

      wlnu = wl*nu
      DO a=1, eNoN
         T1 = N(a) + udNx(a)*tauM
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T2 = wl*(T1*(udNx(b) + (amd + st)*N(b)) 
     2         + tauB*updNx(a)*updNx(b) + nu*NxdNx)
             
            lK(1,a,b) = lK(1,a,b) + wlnu*Nx(1,a)*Nx(1,b) + T2
            lK(2,a,b) = lK(2,a,b) + wlnu*Nx(2,a)*Nx(1,b)
            lK(3,a,b) = lK(3,a,b) + wlnu*Nx(3,a)*Nx(1,b)
            lK(4,a,b) = lK(4,a,b) + wlnu*Nx(1,a)*Nx(2,b)
            lK(5,a,b) = lK(5,a,b) + wlnu*Nx(2,a)*Nx(2,b) + T2
            lK(6,a,b) = lK(6,a,b) + wlnu*Nx(3,a)*Nx(2,b)
            lK(7,a,b) = lK(7,a,b) + wlnu*Nx(1,a)*Nx(3,b)
            lK(8,a,b) = lK(8,a,b) + wlnu*Nx(2,a)*Nx(3,b)
            lK(9,a,b) = lK(9,a,b) + wlnu*Nx(3,a)*Nx(3,b) + T2
         END DO
      END DO
      
      RETURN
      END SUBROUTINE eval3Bbo
!---------------------------------------------------------------------
      PURE SUBROUTINE eval2Bbo(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(bboType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
     
      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,1D0,4D0/)
      INTEGER a, b, Ac
      REAL(KIND=8) wr, T1, T2, amd, wl, udNx(eNoN), res(nsd), tauM, kU, 
     2   nu, tauB, up(nsd), updNx(eNoN), tmpV(nsd), NxdNx, rM(nsd,nsd), 
     3   tauMs, wlnu, ud(nsd), u(nsd), ux(nsd,nsd), bf(nsd), st, 
     4   px(nsd), rho, pRho, pDia
 
      pDia  = eq%mat%D()
      pRho  = eq%mat%rho()
      rho   = eq%fmat%rho()
      nu    = eq%fmat%mu()/rho
      bf(1) = eq%mat%fx() - eq%fmat%fx()
      bf(2) = eq%mat%fy() - eq%fmat%fy()
      st    = 18D0*nu/pDia/pDia/(pRho + 5D-1*rho)
      
      T1    = eq%itg%af*eq%itg%gam*dt
      amd   = eq%itg%am/T1
      wr    = w*(pRho + 5D-1*rho)
      wl    = wr*T1
      
      ud = 0D0
      u  = 0D0
      up = 0D0
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         up(1) = up(1) + N(a)*eq%U%v(1,Ac)
         up(2) = up(2) + N(a)*eq%U%v(2,Ac)

         ud(1) = ud(1) + N(a)*eq%Uf%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%Uf%A%v(2,Ac)
 
         u(1) = u(1) + N(a)*eq%Uf%v(1,Ac)
         u(2) = u(2) + N(a)*eq%Uf%v(2,Ac)
 
         px(1) = px(1) + Nx(1,a)*eq%P%s(Ac)
         px(2) = px(2) + Nx(2,a)*eq%P%s(Ac)

         ux(1,1) = ux(1,1) + Nx(1,a)*eq%Uf%v(1,Ac)
         ux(2,1) = ux(2,1) + Nx(2,a)*eq%Uf%v(1,Ac)
         ux(1,2) = ux(1,2) + Nx(1,a)*eq%Uf%v(2,Ac)
         ux(2,2) = ux(2,2) + Nx(2,a)*eq%Uf%v(2,Ac)
      END DO
 
      bf(1) = st*u(1) + (-rho*px(1) + 5D-1*rho*(ud(1) + up(1)*ux(1,1)
     2   + up(2)*ux(2,1)) + bf(1))/(pRho + 5D-1*rho)
      bf(2) = st*u(2) + (-rho*px(2) + 5D-1*rho*(ud(2) + up(1)*ux(1,2)
     2   + up(2)*ux(2,2)) + bf(2))/(pRho + 5D-1*rho)
       
      u  = up
      ud = 0D0
      ux = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         ud(1) = ud(1) + N(a)*eq%U%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%U%A%v(2,Ac)
 
         ux(1,1) = ux(1,1) + Nx(1,a)*eq%U%v(1,Ac)
         ux(2,1) = ux(2,1) + Nx(2,a)*eq%U%v(1,Ac)
         ux(1,2) = ux(1,2) + Nx(1,a)*eq%U%v(2,Ac)
         ux(2,2) = ux(2,2) + Nx(2,a)*eq%U%v(2,Ac)
      END DO

      kU = u(1)*u(1)*ks(1,1) + u(2)*u(1)*ks(2,1) 
     2   + u(1)*u(2)*ks(1,2) + u(2)*u(2)*ks(2,2)

      tauM = ct(3)/SQRT(ct(1)/dt/dt + ct(2)*kU + st*st)
      tauMs = ct(3)/SQRT(ct(2)*kU + st*st)
      
      res(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + st*u(1) - bf(1)
      res(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + st*u(2) - bf(2)

!     Discontinuity capturing nu (wx , ux)
      tmpV(1) = (u(1)*ux(1,1) + u(2)*ux(1,2))
      tmpV(2) = (u(1)*ux(2,1) + u(2)*ux(2,2))
      
      nu = tmpV(1)*tmpV(1)*ks(1,1) + tmpV(2)*tmpV(1)*ks(2,1) 
     2   + tmpV(1)*tmpV(2)*ks(1,2) + tmpV(2)*tmpV(2)*ks(2,2)
      IF (ISZERO(nu)) nu = eps
      nu = (u(1)*u(1) + u(2)*u(2))/nu
      nu = SQRT(nu*(res(1)*res(1) + res(2)*res(2)))/2D0
      IF (ISZERO(kU)) kU = eps
      T1 = (u(1)*u(1) + u(2)*u(2))/SQRT(kU)/2D0
      IF (nu .GT. T1) nu = T1

!     tauB terms tauB (up.wx , up.ux)
      up = -tauM*res
      tauB = up(1)*up(1)*ks(1,1) + up(2)*up(1)*ks(2,1) 
     2     + up(1)*up(2)*ks(1,2) + up(2)*up(2)*ks(2,2)
      IF (ISZERO(tauB)) tauB = eps
      tauB = ct(3)/SQRT(tauB)
      
      rM(1,1) = nu*(ux(1,1) + ux(1,1)) 
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) 
      rM(1,2) = nu*(ux(1,2) + ux(2,1))
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) 

      tmpV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      tmpV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))
      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a) 
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a)
         T1       = N(a) + udNx(a)*tauM

         lR(1,a) = lR(1,a) + wr*(res(1)*T1 + updNx(a)*tmpV(1) 
     2           + rM(1,1)*Nx(1,a) + rM(1,2)*Nx(2,a))
         lR(2,a) = lR(2,a) + wr*(res(2)*T1 + updNx(a)*tmpV(2)
     2           + rM(2,1)*Nx(1,a) + rM(2,2)*Nx(2,a))
      END DO

      wlnu = wl*nu
      DO a=1, eNoN
         T1 = N(a) + udNx(a)*tauM
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)

            T2 = wl*(T1*(udNx(b) + (amd + st)*N(b)) 
     2         + tauB*updNx(a)*updNx(b) + nu*NxdNx)
             
            lK(1,a,b) = lK(1,a,b) + wlnu*Nx(1,a)*Nx(1,b) + T2
            lK(2,a,b) = lK(2,a,b) + wlnu*Nx(2,a)*Nx(1,b)
            lK(3,a,b) = lK(3,a,b) + wlnu*Nx(1,a)*Nx(2,b)
            lK(4,a,b) = lK(4,a,b) + wlnu*Nx(2,a)*Nx(2,b) + T2
         END DO
      END DO

      RETURN
      END SUBROUTINE eval2Bbo
!---------------------------------------------------------------------
      PURE SUBROUTINE bEvalBbo(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(bboType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)


      RETURN
      END SUBROUTINE bEvalBbo

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_BBOMOD(ns, sct, gg, res)
      TYPE(insType), INTENT(IN) :: ns
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(varType), INTENT(INOUT) :: res

      INTEGER bType
      TYPE(bboType) eq

      cTS   = 0
      eq    = bboType(ns, sct, 'steel', 1)
      bType = IBSET(0,bType_Dir)
      eq%bc(1) = bcType(ns%dmn, "Face_1", gg, bType) ! .,fa_name,.,bType
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
      END SUBROUTINE TEST_BBOMOD
      END MODULE BBOMOD
