!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This is for solving incompressible Navier-Stokes (INS) equations. 
!     For usage and testing see TEST_INSMOD at the end of this file.
!      
!---------------------------------------------------------------------
      MODULE INSMOD
      USE EQMOD
      IMPLICIT NONE

!     Is used to initialize equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 2,
     2   itgO = 1,
     3   lsAl = LS_TYPE_NS,
     5   sym  = "NS")

!     Equation type. Extending solution control 
      TYPE, EXTENDS(eqType) :: insType
!        velocity
         TYPE(gVarType), POINTER :: U => NULL()
!        pressure
         TYPE(gVarType), POINTER :: P => NULL()
!        Two-way coupling force from particle
         TYPE(varType) ::twc
      CONTAINS
!        Setups all structure
         PROCEDURE :: setup => setupIns
!        Evaulate the governing equation
         PROCEDURE :: eval3 => eval3Ins
         PROCEDURE :: eval2 => eval2Ins
!        Evaulate the governing equation boundary contribution
         PROCEDURE :: bEval => bEvalIns
      END TYPE insType

      INTERFACE insType
         PROCEDURE :: newIns, newInsFl
      END INTERFACE insType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newIns(dmn, sct, mName, nBc, ls, itg, cbc) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: mName
      INTEGER, INTENT(IN) :: nBc
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      TYPE(itgType), INTENT(IN), OPTIONAL :: itg
      TYPE(cbcType), INTENT(IN), OPTIONAL :: cbc
      TYPE(insType) :: eq
      
      CALL eq%new(eqSp, dmn, sct, mName, nBc, ls, itg, cbc)

      RETURN
      END FUNCTION newIns
!---------------------------------------------------------------------
      FUNCTION newInsFl(dmn, lst) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(insType) :: eq
      
      CALL eq%new(eqSp, dmn, lst)

      RETURN
      END FUNCTION newInsFl
!---------------------------------------------------------------------
      SUBROUTINE setupIns(eq, var)
      CLASS(insType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
      
      IF (PRESENT(var)) THEN
         IF (SIZE(var) .NE. 2) io%e = "setupIns: Invalid var size"
         IF (var(1)%dof.NE.nsd .OR. var(2)%dof.NE.1) io%e = 
     2      "setupIns: Invalid var number of unknowns"
         eq%var => var
      ELSE
         eq%var(1) = gVarType(nsd,'Velocity',eq%dmn)
         eq%var(2) = gVarType(1,'Pressure',eq%dmn)
      END IF
      eq%U => eq%var(1)
      eq%P => eq%var(2)
      eq%twc = varType(nsd,'TwoWayCoupling',eq%dmn)

      IF (eq%mat%rho() .LE. 0D0)  io%e = "setupIns: rho <= 0"
      IF (eq%mat%mu() .LE. 0D0)   io%e = "setupIns: mu <= 0"
      IF (eq%mat%beta() .GT. 1D0) io%e = "setupIns: beta > 1.0"
      IF (eq%mat%beta() .LT. 0D0) io%e = "setupIns: beta < 0.0"

      RETURN
      END SUBROUTINE setupIns
!---------------------------------------------------------------------
      PURE SUBROUTINE eval3Ins(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(insType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)

      REAL(KIND=8), PARAMETER :: ct(2) = (/1D0,3D0/)
      INTEGER a, b, Ac
      REAL(KIND=8) tauM, tauC, tauB, kssk, kU, nu, rho, T1, T2, T3, 
     2   divU, up(nsd), ua(nsd), udNx(eNoN), updNx(eNoN), uadNx(eNoN), 
     3   rV(nsd), rM(nsd,nsd), NxdNx, amd, wl, wr, wrl, s, ud(nsd), 
     4   u(nsd), p, ux(nsd,nsd), px(nsd), f(nsd)

      rho  = eq%mat%rho()
      nu   = eq%mat%mu()/rho
      f(1) = eq%mat%fx()/rho
      f(2) = eq%mat%fy()/rho
      f(3) = eq%mat%fz()/rho
      s    = eq%mat%A()
     
      wr  = w*rho
      T1  = eq%itg%af*eq%itg%gam*dt
      amd = eq%itg%am/T1
      wrl = wr*T1
      wl  = w*T1
 
      p  = 0D0
      u  = 0D0
      ud = -f
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         Ac = eqN(a)

         p = p + N(a)*eq%P%s(Ac)
               
         px(1) = px(1) + Nx(1,a)*eq%P%s(Ac)
         px(2) = px(2) + Nx(2,a)*eq%P%s(Ac)
         px(3) = px(3) + Nx(3,a)*eq%P%s(Ac)
      
         ud(1) = ud(1) + N(a)*(eq%U%A%v(1,Ac)
     2    + eq%twc%v(1,AC))
         ud(2) = ud(2) + N(a)*(eq%U%A%v(2,Ac)
     2    + eq%twc%v(2,AC))
         ud(3) = ud(3) + N(a)*(eq%U%A%v(3,Ac)
     2    + eq%twc%v(3,AC))

         u(1) = u(1) + N(a)*eq%U%v(1,Ac)
         u(2) = u(2) + N(a)*eq%U%v(2,Ac)
         u(3) = u(3) + N(a)*eq%U%v(3,Ac)
         
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

      IF (ASSOCIATED(eq%dmn%Um)) THEN
         DO a=1, eNoN
            Ac = eqN(a)
            u(1) = u(1) - N(a)*eq%dmn%Um%v(1,Ac)
            u(2) = u(2) - N(a)*eq%dmn%Um%v(2,Ac)
            u(3) = u(3) - N(a)*eq%dmn%Um%v(3,Ac)
         END DO
      END IF

      kU = u(1)*u(1)*ks(1,1) + u(2)*u(1)*ks(2,1) + u(3)*u(1)*ks(3,1) 
     2   + u(1)*u(2)*ks(1,2) + u(2)*u(2)*ks(2,2) + u(3)*u(2)*ks(3,2) 
     3   + u(1)*u(3)*ks(1,3) + u(2)*u(3)*ks(2,3) + u(3)*u(3)*ks(3,3)

      kssk = ks(1,1)*ks(1,1) + ks(2,1)*ks(2,1) + ks(3,1)*ks(3,1) 
     2     + ks(1,2)*ks(1,2) + ks(2,2)*ks(2,2) + ks(3,2)*ks(3,2)
     3     + ks(1,3)*ks(1,3) + ks(2,3)*ks(2,3) + ks(3,3)*ks(3,3)

      tauM = 1D0/SQRT((2D0*ct(1)/dt)**2D0 + kU + ct(2)*nu*nu*kssk + s*s)
      
      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1) 
     2      + u(3)*ux(3,1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2) 
     2      + u(3)*ux(3,2) + s*u(2))
      up(3) = -tauM*(ud(3) + px(3)/rho + u(1)*ux(1,3) + u(2)*ux(2,3) 
     2      + u(3)*ux(3,3) + s*u(3))
            
      tauC = ks(1,1) + ks(2,2) + ks(3,3)
      tauC = 1D0/tauM/tauC/16D0
            
      tauB=up(1)*up(1)*ks(1,1)+ up(2)*up(1)*ks(2,1)+ up(3)*up(1)*ks(3,1)
     2   + up(1)*up(2)*ks(1,2)+ up(2)*up(2)*ks(2,2)+ up(3)*up(2)*ks(3,2)
     4   + up(1)*up(3)*ks(1,3)+ up(2)*up(3)*ks(2,3)+ up(3)*up(3)*ks(3,3)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1D0/SQRT(tauB)

      divU = ux(1,1) + ux(2,2) + ux(3,3) 
   
      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)
      ua(3) = u(3) + up(3)
      
      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
      rV(3) = tauB*(up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))

      rM(1,1) = nu*(ux(1,1) + ux(1,1)) - up(1)*ua(1) + rV(1)*up(1) 
     2        + divU*tauC - p/rho
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) - up(2)*ua(1) + rV(2)*up(1) 
      rM(3,1) = nu*(ux(3,1) + ux(1,3)) - up(3)*ua(1) + rV(3)*up(1) 
      
      rM(1,2) = nu*(ux(1,2) + ux(2,1)) - up(1)*ua(2) + rV(1)*up(2) 
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) - up(2)*ua(2) + rV(2)*up(2) 
     2        + divU*tauC - p/rho
      rM(3,2) = nu*(ux(3,2) + ux(2,3)) - up(3)*ua(2) + rV(3)*up(2) 
      
      rM(1,3) = nu*(ux(1,3) + ux(3,1)) - up(1)*ua(3) + rV(1)*up(3) 
      rM(2,3) = nu*(ux(2,3) + ux(3,2)) - up(2)*ua(3) + rV(2)*up(3) 
      rM(3,3) = nu*(ux(3,3) + ux(3,3)) - up(3)*ua(3) + rV(3)*up(3) 
     2        + divU*tauC - p/rho

      rV(1) = ud(1) + ua(1)*(s+ux(1,1)) + ua(2)*ux(2,1) + ua(3)*ux(3,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*(s+ux(2,2)) + ua(3)*ux(3,2)
      rV(3) = ud(3) + ua(1)*ux(1,3) + ua(2)*ux(2,3) + ua(3)*(s+ux(3,3))
      
      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)  + u(3)*Nx(3,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a) + up(3)*Nx(3,a)
         uadNx(a) = updNx(a) + udNx(a)

         lR(1,a) = lR(1,a) + wr*(rV(1)*N(a) + rM(1,1)*Nx(1,a) 
     2      + rM(1,2)*Nx(2,a) + rM(1,3)*Nx(3,a))

         lR(2,a) = lR(2,a) + wr*(rV(2)*N(a) + rM(2,1)*Nx(1,a) 
     2      + rM(2,2)*Nx(2,a) + rM(2,3)*Nx(3,a))
         
         lR(3,a) = lR(3,a) + wr*(rV(3)*N(a) + rM(3,1)*Nx(1,a) 
     2      + rM(3,2)*Nx(2,a) + rM(3,3)*Nx(3,a))
         
         lR(4,a) = lR(4,a) + w*(N(a)*divU - updNx(a))
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            rM(1,1) = Nx(1,a)*Nx(1,b)
            rM(2,1) = Nx(2,a)*Nx(1,b)
            rM(3,1) = Nx(3,a)*Nx(1,b)
            rM(1,2) = Nx(1,a)*Nx(2,b)
            rM(2,2) = Nx(2,a)*Nx(2,b)
            rM(3,2) = Nx(3,a)*Nx(2,b)
            rM(1,3) = Nx(1,a)*Nx(3,b)
            rM(2,3) = Nx(2,a)*Nx(3,b)
            rM(3,3) = Nx(3,a)*Nx(3,b)

            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)
 
            T1 = nu*NxdNx + tauB*updNx(a)*updNx(b)
     2         + N(a)*((s+amd)*N(b) + uadNx(b))
     3         + tauM*uadNx(a)*(udNx(b) + (s+amd)*N(b))

            T2 = tauM*udNx(a)
            T3 = tauM*(amd*N(b) + udNx(b))
!           dM/dU
            lK(1,a,b)  = lK(1,a,b)  + wrl*((nu + tauC)*rM(1,1) + T1)
            lK(2,a,b)  = lK(2,a,b)  + wrl*(nu*rM(2,1) + tauC*rM(1,2))
            lK(3,a,b)  = lK(3,a,b)  + wrl*(nu*rM(3,1) + tauC*rM(1,3))
            
            lK(5,a,b)  = lK(5,a,b)  + wrl*(nu*rM(1,2) + tauC*rM(2,1))
            lK(6,a,b)  = lK(6,a,b)  + wrl*((nu + tauC)*rM(2,2) + T1)
            lK(7,a,b)  = lK(7,a,b)  + wrl*(nu*rM(3,2) + tauC*rM(2,3))
            
            lK(9,a,b)  = lK(9,a,b)  + wrl*(nu*rM(1,3) + tauC*rM(3,1))
            lK(10,a,b) = lK(10,a,b) + wrl*(nu*rM(2,3) + tauC*rM(3,2))
            lK(11,a,b) = lK(11,a,b) + wrl*((nu + tauC)*rM(3,3) + T1)
!           dM/dP 
            lK(4,a,b)  = lK(4,a,b)  - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2)
            lK(8,a,b)  = lK(8,a,b)  - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2)
            lK(12,a,b) = lK(12,a,b) - wl*(Nx(3,a)*N(b) - Nx(3,b)*T2)
!           dC/dU
            lK(13,a,b) = lK(13,a,b) + wl*(Nx(1,b)*N(a) + Nx(1,a)*T3)
            lK(14,a,b) = lK(14,a,b) + wl*(Nx(2,b)*N(a) + Nx(2,a)*T3)
            lK(15,a,b) = lK(15,a,b) + wl*(Nx(3,b)*N(a) + Nx(3,a)*T3)
!           dC/dP
            lK(16,a,b) = lK(16,a,b) + wl*(tauM*NxdNx)/rho
         END DO
      END DO
      
      RETURN
      END SUBROUTINE eval3Ins
!---------------------------------------------------------------------
      PURE SUBROUTINE eval2Ins(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(insType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
     
      REAL(KIND=8), PARAMETER :: ct(2) = (/1D0,3D0/)
      INTEGER a, b, Ac
      REAL(KIND=8) tauM, tauC, tauB, kssk, kU, nu, rho, T1, T2, T3, 
     2   divU, up(nsd), ua(nsd), udNx(eNoN), updNx(eNoN), uadNx(eNoN), 
     3   rV(nsd), rM(nsd,nsd), NxdNx, amd, wl, wr, wrl, s, ud(nsd), 
     4   u(nsd), p, ux(nsd,nsd), px(nsd), f(nsd)

      rho  = eq%mat%rho()
      nu   = eq%mat%mu()/rho
      f(1) = eq%mat%fx()/rho
      f(2) = eq%mat%fy()/rho
      s    = eq%mat%A()
      
      wr  = w*rho
      T1  = eq%itg%af*eq%itg%gam*dt
      amd = eq%itg%am/T1
      wrl = wr*T1
      wl  = w*T1
 
      p  = 0D0
      u  = 0D0
      ud = -f
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         Ac = eqN(a)

         p = p + N(a)*eq%P%s(Ac)
               
         px(1) = px(1) + Nx(1,a)*eq%P%s(Ac)
         px(2) = px(2) + Nx(2,a)*eq%P%s(Ac)
      
         ud(1) = ud(1) + N(a)*eq%U%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%U%A%v(2,Ac)

         u(1) = u(1) + N(a)*eq%U%v(1,Ac)
         u(2) = u(2) + N(a)*eq%U%v(2,Ac)
         
         ux(1,1) = ux(1,1) + Nx(1,a)*eq%U%v(1,Ac)
         ux(2,1) = ux(2,1) + Nx(2,a)*eq%U%v(1,Ac)
         ux(1,2) = ux(1,2) + Nx(1,a)*eq%U%v(2,Ac)
         ux(2,2) = ux(2,2) + Nx(2,a)*eq%U%v(2,Ac)
      END DO

      IF (ASSOCIATED(eq%dmn%Um)) THEN
         DO a=1, eNoN
            Ac = eqN(a)
            u(1) = u(1) - N(a)*eq%dmn%Um%v(1,Ac)
            u(2) = u(2) - N(a)*eq%dmn%Um%v(2,Ac)
         END DO
      END IF

      kU = u(1)*u(1)*ks(1,1) + u(2)*u(1)*ks(2,1)
     2   + u(1)*u(2)*ks(1,2) + u(2)*u(2)*ks(2,2)

      kSsk = ks(1,1)*ks(1,1) + ks(2,1)*ks(2,1)
     2     + ks(1,2)*ks(1,2) + ks(2,2)*ks(2,2)
         
      tauM = 1D0/SQRT((2D0*ct(1)/dt)**2D0 + kU + ct(2)*nu*nu*kssk + s*s)
            
      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + s*u(2))
            
      tauC = ks(1,1) + ks(2,2)
      tauC = 1D0/tauM/tauC/16D0
            
      tauB = up(1)*up(1)*ks(1,1) + up(2)*up(1)*ks(2,1) 
     2     + up(1)*up(2)*ks(1,2) + up(2)*up(2)*ks(2,2)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1D0/SQRT(tauB)

      divU = ux(1,1) + ux(2,2)
   
      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)
      
      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))

      rM(1,1) = nu*(ux(1,1) + ux(1,1)) - up(1)*ua(1) + rV(1)*up(1) 
     2        + divU*tauC - p/rho
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) - up(2)*ua(1) + rV(2)*up(1) 
      
      rM(1,2) = nu*(ux(1,2) + ux(2,1)) - up(1)*ua(2) + rV(1)*up(2) 
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) - up(2)*ua(2) + rV(2)*up(2) 
     2        + divU*tauC - p/rho
      
      rV(1) = ud(1) + ua(1)*(s + ux(1,1)) + ua(2)*ux(2,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*(s + ux(2,2))
 
      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a)
         uadNx(a) = updNx(a) + udNx(a)

         lR(1,a) = lR(1,a) + wr*(rV(1)*N(a) + rM(1,1)*Nx(1,a) 
     2      + rM(1,2)*Nx(2,a))

         lR(2,a) = lR(2,a) + wr*(rV(2)*N(a) + rM(2,1)*Nx(1,a) 
     2      + rM(2,2)*Nx(2,a))
         
         lR(3,a) = lR(3,a) + w*(N(a)*divU - updNx(a))
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            rM(1,1) = Nx(1,a)*Nx(1,b)
            rM(2,1) = Nx(2,a)*Nx(1,b)
            rM(1,2) = Nx(1,a)*Nx(2,b)
            rM(2,2) = Nx(2,a)*Nx(2,b)

            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)
        
            T1 = nu*NxdNx + tauB*updNx(a)*updNx(b)
     2         + N(a)*((s+amd)*N(b) + uadNx(b))
     3         + tauM*uadNx(a)*(udNx(b) + (s+amd)*N(b))

            T2 = tauM*udNx(a)
            T3 = tauM*(amd*N(b) + udNx(b))
!           dM/dU
            lK(1,a,b)  = lK(1,a,b) + wrl*((nu + tauC)*rM(1,1) + T1)
            lK(2,a,b)  = lK(2,a,b) + wrl*(nu*rM(2,1) + tauC*rM(1,2))
            
            lK(4,a,b)  = lK(4,a,b) + wrl*(nu*rM(1,2) + tauC*rM(2,1))
            lK(5,a,b)  = lK(5,a,b) + wrl*((nu + tauC)*rM(2,2) + T1)
!           dM/dP
            lK(3,a,b)  = lK(3,a,b) - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2)
            lK(6,a,b)  = lK(6,a,b) - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2)
!           dC/dU
            lK(7,a,b) = lK(7,a,b)  + wl*(Nx(1,b)*N(a) + Nx(1,a)*T3)
            lK(8,a,b) = lK(8,a,b)  + wl*(Nx(2,b)*N(a) + Nx(2,a)*T3)
!           dC/dP
            lK(9,a,b) = lK(9,a,b)  + wl*tauM*NxdNx/rho
         END DO
      END DO

      RETURN
      END SUBROUTINE eval2Ins
!---------------------------------------------------------------------
      PURE SUBROUTINE bEvalIns(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(insType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)

      INTEGER a, b, Ac
      REAL(KIND=8) T1, wl, hc(nsd), udn, u(nsd), beta, rho

      rho  = eq%mat%rho()
      beta = eq%mat%beta()
      wl   = w*eq%itg%af*eq%itg%gam*dt

      u = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         u  = u + N(a)*eq%U%v(:,Ac)
      END DO

      IF (ASSOCIATED(eq%dmn%Um)) THEN
         DO a=1, eNoN
            Ac = eqN(a)
            u = u - N(a)*eq%dmn%Um%v(:,Ac)
         END DO
      END IF

      udn = SUM(u*nV)
      udn = rho*beta*(udn - ABS(udn))/2D0
      hc  = h*nV + udn*u
 
!     Here the loop is started for constructing left and right hand side
      IF (nsd .EQ. 2) THEN
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - w*N(a)*hc(1)
            lR(2,a) = lR(2,a) - w*N(a)*hc(2)
            DO b=1, eNoN
               T1        = wl*N(a)*N(b)*udn
               lK(1,a,b) = lK(1,a,b) - T1
               lK(5,a,b) = lK(5,a,b) - T1
            END DO
         END DO
      ELSE
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - w*N(a)*hc(1)
            lR(2,a) = lR(2,a) - w*N(a)*hc(2)
            lR(3,a) = lR(3,a) - w*N(a)*hc(3)
            DO b=1, eNoN
               T1 = wl*N(a)*N(b)*udn
               lK(1,a,b)  = lK(1,a,b)  - T1
               lK(6,a,b)  = lK(6,a,b)  - T1
               lK(11,a,b) = lK(11,a,b) - T1
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE bEvalIns

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_INSMOD(dmn, sct, gg, eq, res)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(insType), INTENT(OUT) :: eq
      TYPE(varType), INTENT(INOUT) :: res

      INTEGER bType

      cTS   = 0
      eq    = insType(dmn, sct, 'water', 1) ! ., ., mat_name, nBc
      bType = IBSET(0,bType_Dir)
      eq%bc(1) = bcType(dmn, "face_1", gg, bType) ! ., fa_name, ., bType
      CALL eq%ini()
      DO cTS=1, nTS
         DO WHILE(.NOT.eq%satisfied())
            CALL eq%solve()
         END DO
         CALL eq%update()
      END DO
      res%s = eq%U%v(2,:)

      RETURN
      END SUBROUTINE TEST_INSMOD
!---------------------------------------------------------------------
!     To test coupled BC.
      SUBROUTINE TEST_CBCMOD(dmn, sct, gg, res)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(varType), INTENT(INOUT) :: res

      LOGICAL fl
      INTEGER bType
      TYPE(cbcType) cbc
      TYPE(fileType) f
      TYPE(insType) eq

      INQUIRE(FILE='cplBC.exe', EXIST=fl)
      IF (.NOT.fl) io%e = "cplBC.exe is needed in the current directory"
      
      f = fileType('0D.ini','ascii')
      CALL f%open('w')
      WRITE(f%id(),*) 1D1
      CALL f%close()

      cbc = cbcType('semi-implicit','./cplBC.exe', f, 1) ! .,.,.,nX
      io%o = "newCbc: "//CLR("(PASSED)",3)
      CALL cbc%addBc('Face_1', cbc_Neu)
      io%o = "cbc%addBc: "//CLR("(PASSED)",3)

      cTS   = 0
      eq    = insType(dmn, sct, 'water', 1, cbc=cbc) !.,.,mat_name,nBc,.
      bType = IBSET(0,bType_Neu)
      bType = IBSET(bType,bType_cpl)
      eq%bc(1) = bcType(dmn, "face_1", gg, bType, cbc=cbc)
      CALL eq%ini()
      DO cTS=1, nTS
         CALL eq%update()
         DO WHILE(.NOT.eq%satisfied())
            CALL eq%solve()
         END DO
         CALL eq%save()
      END DO
      res%s = eq%U%v(2,:)
      CALL cbc%free()
      CALL eq%free()

      RETURN
      END SUBROUTINE TEST_CBCMOD
      END MODULE INSMOD
