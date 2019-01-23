!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This is for solving linear elastodynamics (led) equations. 
!     For usage and testing see TEST_LEDMOD at the end of this file.
!      
!---------------------------------------------------------------------
      MODULE LEDMOD
      USE EQMOD
      IMPLICIT NONE

!     Is used to initialize equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 1,
     2   itgO = 2,
     3   lsAl = LS_TYPE_CG,
     4   sym  = "LE")

!     Equation type. Extending solution control 
      TYPE, EXTENDS(eqType) :: ledType
         PRIVATE
!        Whether this is solving for mesh motion
         LOGICAL :: isMsh = .FALSE. 
!        Displacement, velocity, ...
         TYPE(gVarType), POINTER :: U => NULL()
      CONTAINS
!        Setups all structure
         PROCEDURE :: setup => setupLed
!        Evaulate the governing equation
         PROCEDURE :: eval3 => eval3Led
         PROCEDURE :: eval2 => eval2Led
!        Overridden procedures (no implementation)
         PROCEDURE :: bEval => bEvalLed
      END TYPE ledType

      INTERFACE ledType
         PROCEDURE :: newLed, newLedFl
      END INTERFACE ledType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newLed(dmn, sct, mName, nBc, ls, isMsh, g, flag) 
     2   RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: mName
      INTEGER, INTENT(IN) :: nBc
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      LOGICAL, INTENT(IN), OPTIONAL :: isMsh
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: g
      LOGICAL, INTENT(IN), OPTIONAL :: flag(dmn%nMsh)
      TYPE(ledType) :: eq
 
      CALL eq%new(eqSp, dmn, sct, mName, nBc, ls)

      IF (PRESENT(isMsh)) eq%isMsh = isMsh
      IF (PRESENT(g).NEQV.PRESENT(flag)) io%e = "newLed: g and flag"//
     2   " must be present at the same time"
      IF (PRESENT(g)) THEN
         IF (.NOT.PRESENT(isMsh)) io%e = "newLedFl: g/flag can be"//
     2      " defined only when isMsh=.TRUE."
         eq%vc = vcType(g, flag)
      END IF

      RETURN
      END FUNCTION newLed
!---------------------------------------------------------------------
      FUNCTION newLedFl(dmn, lst, isMsh, g, flag) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst
      LOGICAL, INTENT(IN), OPTIONAL :: isMsh
      TYPE(gVarType), INTENT(IN), OPTIONAL :: g
      LOGICAL, INTENT(IN), OPTIONAL :: flag(dmn%nMsh)
      TYPE(ledType) :: eq

      CALL eq%new(eqSp, dmn, lst)

      IF (PRESENT(isMsh)) eq%isMsh = isMsh
      IF (PRESENT(g).NEQV.PRESENT(flag)) io%e = "newLedFl: g and flag"//
     2   " must be present at the same time"
      IF (PRESENT(g)) THEN
         IF (.NOT.PRESENT(isMsh)) io%e = "newLedFl: g/flag can be"//
     2      " defined only when isMsh=.TRUE."
         eq%vc = vcType(g, flag)
      END IF

      RETURN
      END FUNCTION newLedFl
!---------------------------------------------------------------------
      SUBROUTINE setupLed(eq, var)
      CLASS(ledType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
      
      IF (PRESENT(var)) THEN
         IF (SIZE(var) .NE. 1) io%e = "setupLed: Invalid var size"
         IF (var(1)%dof .NE. nsd) io%e = "setupLed: Invalid var nU"
         eq%var => var
      ELSE
         eq%var(1) = gVarType(nsd,'Velocity',eq%dmn)
      END IF
      eq%U => eq%var(1)
      IF (eq%isMsh) THEN
         IF (ASSOCIATED(eq%dmn%Um)) io%e = "MSH equation can be"//
     2      " specified only once"
         eq%dmn%Um => eq%var(1)
         eq%conf   = eqConf_old
      END IF
     
      IF (eq%mat%rho() .LE. 0D0) io%e = "setupLed: rho <= 0"
      IF (eq%mat%E()   .LE. 0D0) io%e = "setupLed: E <= 0"
      IF (eq%mat%nu() .GE. 5D-1) io%e = "setupLed: nu >= 0.5"
      IF (eq%mat%nu()  .LE. 0D0) io%e = "setupLed: nu <= 0.0"

      RETURN
      END SUBROUTINE setupLed
!---------------------------------------------------------------------
      PURE SUBROUTINE eval3Led(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(ledType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
      
      INTEGER a, b, Ac
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd, wl, 
     2   lDm, ed(6), ud(nsd), f(nsd), wr
 
      wr = w
      IF (eq%isMsh) wr = w/J

      rho  = eq%mat%rho()
      elM  = eq%mat%E()
      nu   = eq%mat%nu() 
      f(1) = eq%mat%fx()
      f(2) = eq%mat%fy()
      f(3) = eq%mat%fz()
      
      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      lDm    = lambda/mu
      T1     = eq%itg%af*eq%itg%beta*dt*dt
      amd    = eq%itg%am/T1*rho
      wl     = wr*T1*mu

      ed = 0D0
      ud = -f
      DO a=1, eNoN
         Ac = eqN(a)
         ud(1) = ud(1) + N(a)*eq%U%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%U%A%v(2,Ac)
         ud(3) = ud(3) + N(a)*eq%U%A%v(3,Ac)

         ed(1) = ed(1) + Nx(1,a)*eq%U%D%v(1,Ac)
         ed(2) = ed(2) + Nx(2,a)*eq%U%D%v(2,Ac)
         ed(3) = ed(3) + Nx(3,a)*eq%U%D%v(3,Ac)
         ed(4) = ed(4) + Nx(3,a)*eq%U%D%v(2,Ac) + Nx(2,a)*eq%U%D%v(3,Ac)
         ed(5) = ed(5) + Nx(3,a)*eq%U%D%v(1,Ac) + Nx(1,a)*eq%U%D%v(3,Ac)
         ed(6) = ed(6) + Nx(2,a)*eq%U%D%v(1,Ac) + Nx(1,a)*eq%U%D%v(2,Ac)
      END DO
      IF (eq%isMsh) THEN
         DO a=1, eNoN
            Ac = eqN(a)
            ed(1) = ed(1) - Nx(1,a)*eq%U%Do%v(1,Ac)
            ed(2) = ed(2) - Nx(2,a)*eq%U%Do%v(2,Ac)
            ed(3) = ed(3) - Nx(3,a)*eq%U%Do%v(3,Ac)
            ed(4) = ed(4) - Nx(3,a)*eq%U%Do%v(2,Ac) 
     2                    - Nx(2,a)*eq%U%Do%v(3,Ac)
            ed(5) = ed(5) - Nx(3,a)*eq%U%Do%v(1,Ac) 
     2                    - Nx(1,a)*eq%U%Do%v(3,Ac)
            ed(6) = ed(6) - Nx(2,a)*eq%U%Do%v(1,Ac) 
     2                    - Nx(1,a)*eq%U%Do%v(2,Ac)
         END DO
      END IF
      divD = lambda*(ed(1) + ed(2) + ed(3))
 
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + wr*(rho*N(a)*ud(1) + mu*Nx(2,a)*ed(6) 
     2      + mu*Nx(3,a)*ed(5) + Nx(1,a)*(2D0*mu*ed(1) + divD))

         lR(2,a) = lR(2,a) + wr*(rho*N(a)*ud(2) + mu*Nx(1,a)*ed(6) 
     2      + mu*Nx(3,a)*ed(4) + Nx(2,a)*(2D0*mu*ed(2) + divD))
         
         lR(3,a) = lR(3,a) + wr*(rho*N(a)*ud(3) + mu*Nx(1,a)*ed(5) 
     2      + mu*Nx(2,a)*ed(4) + Nx(3,a)*(2D0*mu*ed(3) + divD))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T1 = amd*N(a)*N(b)/mu + NxdNx

            lK(1,a,b) = lK(1,a,b) + wl*(T1 + (1D0+lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b) 
     2         + Nx(2,a)*Nx(1,b))

            lK(3,a,b) = lK(3,a,b) + wl*(lDm*Nx(1,a)*Nx(3,b) 
     2         + Nx(3,a)*Nx(1,b))

            lK(4,a,b) = lK(4,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b) 
     2         + Nx(1,a)*Nx(2,b))

            lK(5,a,b) = lK(5,a,b) + wl*(T1 + (1D0+lDm)*Nx(2,a)*Nx(2,b))
            
            lK(6,a,b) = lK(6,a,b) + wl*(lDm*Nx(2,a)*Nx(3,b) 
     2         + Nx(3,a)*Nx(2,b))

            lK(7,a,b) = lK(7,a,b) + wl*(lDm*Nx(3,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(3,b))

            lK(8,a,b) = lK(8,a,b) + wl*(lDm*Nx(3,a)*Nx(2,b) 
     2         + Nx(2,a)*Nx(3,b))
            
            lK(9,a,b) = lK(9,a,b) + wl*(T1 + (1D0+lDm)*Nx(3,a)*Nx(3,b))
         END DO
      END DO
      
      RETURN
      END SUBROUTINE eval3Led
!---------------------------------------------------------------------
      PURE SUBROUTINE eval2Led(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(ledType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
      
      INTEGER a, b, Ac
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd, wl, 
     2   lDm, ed(3), ud(nsd), f(nsd), wr

      wr = w
      IF (eq%isMsh) wr = w/J

      rho  = eq%mat%rho()
      elM  = eq%mat%E()
      nu   = eq%mat%nu() 
      f(1) = eq%mat%fx()
      f(2) = eq%mat%fy()
      
      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      lDm    = lambda/mu
      T1     = eq%itg%af*eq%itg%beta*dt*dt
      amd    = eq%itg%am/T1*rho
      wl     = wr*T1*mu

      ed = 0D0
      ud = -f
      DO a=1, eNoN
         Ac = eqN(a)
         ud(1) = ud(1) + N(a)*eq%U%A%v(1,Ac)
         ud(2) = ud(2) + N(a)*eq%U%A%v(2,Ac)

         ed(1) = ed(1) + Nx(1,a)*eq%U%D%v(1,Ac)
         ed(2) = ed(2) + Nx(2,a)*eq%U%D%v(2,Ac)
         ed(3) = ed(3) + Nx(2,a)*eq%U%D%v(1,Ac) + Nx(1,a)*eq%U%D%v(2,Ac)
      END DO
      IF (eq%isMsh) THEN
         DO a=1, eNoN
            Ac = eqN(a)
            ed(1) = ed(1) - Nx(1,a)*eq%U%Do%v(1,Ac)
            ed(2) = ed(2) - Nx(2,a)*eq%U%Do%v(2,Ac)
            ed(3) = ed(3) - Nx(2,a)*eq%U%Do%v(1,Ac) 
     2                    - Nx(1,a)*eq%U%Do%v(2,Ac)
         END DO
      END IF
      divD = lambda*(ed(1) + ed(2))
 
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + wr*(rho*N(a)*ud(1) + mu*Nx(2,a)*ed(3) 
     2                     + Nx(1,a)*(2D0*mu*ed(1) + divD))

         lR(2,a) = lR(2,a) + wr*(rho*N(a)*ud(2) + mu*Nx(1,a)*ed(3) 
     2                     + Nx(2,a)*(2D0*mu*ed(2) + divD))
         
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)
            
            T1 = amd*N(a)*N(b)/mu + NxdNx
            
            lK(1,a,b) = lK(1,a,b) + wl*(T1 + (1D0+lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b) 
     2         + Nx(2,a)*Nx(1,b))

            lK(3,a,b) = lK(3,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b) 
     2         + Nx(1,a)*Nx(2,b))

            lK(4,a,b) = lK(3,a,b) + wl*(T1 + (1D0+lDm)*Nx(2,a)*Nx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE eval2Led
!---------------------------------------------------------------------
      PURE SUBROUTINE bEvalLed(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(ledType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)


      RETURN
      END SUBROUTINE bEvalLed

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_LEDMOD(dmn, sct, gg, res)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(varType), INTENT(INOUT) :: res

      INTEGER bType
      TYPE(ledType) eq

      cTS   = 0
      eq    = ledType(dmn, sct, 'Steel', 1) ! ., ., mat_name, nBc
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
      END SUBROUTINE TEST_LEDMOD
      END MODULE LEDMOD
