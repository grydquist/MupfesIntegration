!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This is for solving incompressible Navier-Stokes (FSI) equations. 
!     For usage and testing see TEST_FSIMOD at the end of this file.
!      
!---------------------------------------------------------------------
      MODULE FSIMOD
      USE INSMOD
      USE SVKMOD
      USE LEDMOD
      IMPLICIT NONE

!     Is used to initialize equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 2,
     2   itgO = 1,
     3   lsAl = LS_TYPE_GMRES,
     4   sym  = "FS")

!     Equation type. Extending solution control 
      TYPE, EXTENDS(eqType) :: fsiType
         PRIVATE
!        velocity
         TYPE(gVarType), POINTER :: U => NULL()
!        pressure
         TYPE(gVarType), POINTER :: P => NULL()
!        The fluid equation
         TYPE(insType) ins
!        The solid equation
         TYPE(svkType) svk
      CONTAINS
!        Setups all structure
         PROCEDURE :: setup => setupFsi
!        Overridden procedures (no implementation)
         PROCEDURE :: eval2 => eval2Fsi
         PROCEDURE :: eval3 => eval3Fsi
         PROCEDURE :: bEval => bEvalFsi
      END TYPE fsiType

      INTERFACE fsiType
         PROCEDURE :: newFsi, newFsiFl
      END INTERFACE fsiType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newFsi(dmn, sct, fmName, smName, nBc, isFluid, ls, 
     2   cbc) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: fmName, smName
      INTEGER, INTENT(IN) :: nBc
      LOGICAL, INTENT(IN) :: isFluid(dmn%nMsh)
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      TYPE(cbcType), INTENT(IN), OPTIONAL :: cbc
      TYPE(fsiType) :: eq
      
      TYPE(lsType) lsTmp
      
      CALL eq%new(eqSp, dmn, sct, fmName, nBc, ls=ls, cbc=cbc)
 
      eq%nSub = 2
      ALLOCATE(eq%sub(eq%nSub), eq%subPtr(dmn%nMsh))
     
      lsTmp%nU  = nsd + 1
      ALLOCATE(eq%sub(1)%s, SOURCE=insType(dmn, eq%sctType, fmName, 0, 
     2   lsTmp, eq%itg))
      
      lsTmp%nU  = nsd
      ALLOCATE(eq%sub(2)%s, SOURCE=svkType(dmn, eq%sctType, smName, 0, 
     2   lsTmp, eq%itg))

      WHERE(isFluid)
         eq%subPtr = 1
      ELSEWHERE
         eq%subPtr = 2
      END WHERE

      RETURN
      END FUNCTION newFsi
!---------------------------------------------------------------------
      FUNCTION newFsiFl(dmn, lst) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(fsiType) :: eq
      
      INTEGER i, n, iM
      CHARACTER(LEN=stdL) ctmp
      TYPE(lstType), POINTER :: lPtr
      TYPE(lsType) ls
      
      CALL eq%new(eqSp, dmn, lst)
      
      eq%nSub = 2
      ALLOCATE(eq%sub(eq%nSub), eq%subPtr(dmn%nMsh))
      lPtr  => lst%get(ctmp,"Material",1)
      ls%nU  = nsd + 1
      ALLOCATE(eq%sub(1)%s, SOURCE=insType(dmn, eq%sctType, ctmp, 0, ls,
     2   eq%itg))
      
      lPtr  => lst%get(ctmp,"Solid material",1)
      ls%nU  = nsd
      ALLOCATE(eq%sub(2)%s, SOURCE=svkType(dmn, eq%sctType, ctmp, 0, ls,
     2   eq%itg))

      eq%subPtr = 1
      n = lst%srch("Solid domain",1)
      DO i=1, n
         lPtr  => lst%get(ctmp,"Solid domain",i)
         DO iM=1, dmn%nMsh
            IF (dmn%msh(iM)%name .EQ. ctmp) EXIT
         END DO
         IF (iM .GT. dmn%nMsh) io%e = "newFsiFl: Unable to find mesh <"
     2      //TRIM(ctmp)//">"
         eq%subPtr(iM) = 2
      END DO

      RETURN
      END FUNCTION newFsiFl
!---------------------------------------------------------------------
      SUBROUTINE setupFsi(eq, var)
      CLASS(fsiType), INTENT(INOUT), TARGET :: eq
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
      eq%conf = eqConf_NA
     
      CALL eq%sub(1)%s%setup(eq%var)
      CALL eq%sub(2)%s%setup(eq%var(1:1))

      RETURN
      END SUBROUTINE setupFsi
!---------------------------------------------------------------------
      SUBROUTINE eval3Fsi(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(fsiType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)

      io%e = "eval3Fsi: this routine is not supposed to be called" 
      
      RETURN
      END SUBROUTINE eval3Fsi
!---------------------------------------------------------------------
      SUBROUTINE eval2Fsi(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(fsiType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
 
      io%e = "eval2Fsi: this routine is not supposed to be called" 

      RETURN
      END SUBROUTINE eval2Fsi
!---------------------------------------------------------------------
      SUBROUTINE bEvalFsi(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(fsiType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)

      io%e = "eval2Fsi: this routine is not supposed to be called" 

      RETURN
      END SUBROUTINE bEvalFsi

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_FSIMOD(dmn, sct, gg, res)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(varType), INTENT(INOUT) :: res

      INTEGER bType
      TYPE(fsiType) eq
      TYPE(ledType) mm
      LOGICAL isFluid(dmn%nMsh), flag(dmn%nMsh)

      IF (dmn%nMsh .NE. 2) THEN
         io%w = "This test case only works with -ts2"
         RETURN
      END IF

      dt       = dt/1D1
      isFluid  = (/.TRUE.,.FALSE./)
      flag     = (/.FALSE.,.TRUE./)
      cTS      = 0
      eq       = fsiType(dmn, sct, 'water', 'steel', 1, isFluid)
      bType    = IBSET(0,bType_Dir)
      eq%bc(1) = bcType(dmn, "face_1", gg, bType) ! ., fa_name, ., bType
      CALL eq%ini()
      mm = ledType(dmn, sct, 'steel', 0, isMsh=.TRUE., g=eq%var(1), 
     2   flag=flag) ! .,.,.,.,.,FSI_var, wh-msh-to-impose
      CALL mm%ini()
      DO cTS=1, nTS
         CALL eq%update()
         CALL mm%update()
         DO WHILE(.NOT.eq%satisfied())
            CALL eq%solve()
         END DO
         DO WHILE(.NOT.mm%satisfied())
            CALL mm%solve()
         END DO
      END DO
      dt    = dt*1D1
      res%s = eq%U%v(2,:)
      CALL eq%free()

      RETURN
      END SUBROUTINE TEST_FSIMOD
      END MODULE FSIMOD
