!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This class is multidomain simulation and communication with 0D
!     solver cplBC. For usage and testing see TEST_CBCMOD at the end of
!     INS.f file. 
!      
!---------------------------------------------------------------------
      MODULE CBCMOD
      USE ITGMOD
      IMPLICIT NONE
      INCLUDE "cplBC.h"

!     Differenty type of coupling for cbc
!     Implicit, semi-implicit, and explicit
      INTEGER, PARAMETER, PRIVATE :: cbc_I = 401, cbc_SI = 402, 
     2   cbc_E = 403
!     File name for communication between 0D and 3D
      CHARACTER(LEN=*), PARAMETER, PRIVATE :: cbc_cmName = 
     2   ".CPLBC_0D_3D.tmp"

!     For coupled 0D-3D problems
      TYPE cbcType
         PRIVATE
!        Weather the calcDer is called before
         LOGICAL :: called = .FALSE.
!        Weather it is initialized
         LOGICAL :: izd = .FALSE.
!        Number of coupled faces
         INTEGER :: nFa = 0
!        Number of unknowns in the 0D domain
         INTEGER :: nX = 0
!        Implicit/Explicit/Semi-implicit schemes
         INTEGER :: schm
!        Path to the 0D code binary file 
         CHARACTER(LEN=stdL) :: binPath
!        New time step unknowns in the 0D domain
         REAL(KIND=8), ALLOCATABLE, PUBLIC :: xn(:)
!        Old time step unknowns in the 0D domain
         REAL(KIND=8), ALLOCATABLE :: xo(:)
!        Data structure used for communicating with 0D code
         TYPE(cplFaceType), ALLOCATABLE :: fa(:)
      CONTAINS
!        Whether the class is initialized
         PROCEDURE :: isIzd => isIzdCbc
!        Whether the scheme is explicit
         PROCEDURE :: isE => isECbc
!        Adds a new coupled BC
         PROCEDURE :: addBc => addBcCbc
!        Solves 0D domain
         PROCEDURE :: solve => solveCbc
!        Returns 0D solution 
         PROCEDURE :: get => getCbc
!        To deallocate the structure
         PROCEDURE :: free => freeCbc
!        Calls Lower Dimensional model
         PROCEDURE, PRIVATE :: callLD => callLDCbc
!        Calculates derivative R = dP/dQ for Neu BCs
         PROCEDURE, PRIVATE :: calcDer => calcDerCbc
      END TYPE cbcType

      INTERFACE cbcType
         PROCEDURE :: newCbc, newCbcFl
      END INTERFACE cbcType

      CONTAINS 
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newCbc(schm, binPath, iniF, nX) RESULT(cbc)
      CHARACTER(LEN=*), INTENT(IN) :: schm, binPath
      TYPE(fileType), INTENT(INOUT) :: iniF
      INTEGER, INTENT(IN) :: nX
      TYPE(cbcType) :: cbc

      SELECT CASE(schm)
      CASE('implicit')
         cbc%schm = cbc_I
      CASE('semi-implicit')
         cbc%schm = cbc_SI
      CASE('explicit')
         cbc%schm = cbc_E
      CASE DEFAULT
         io%e = "Undefined cbc%schm: "//schm
      END SELECT

      cbc%izd = .TRUE.
      cbc%nX  = nX
      ALLOCATE(cbc%xo(cbc%nX), cbc%xn(cbc%nX))
      cbc%xo = 0D0

      IF (cm%mas()) THEN
         CALL iniF%open('r')
         CALL iniF%rw(cbc%xo)
         CALL iniF%close()
      END IF
      CALL cm%bcast(cbc%xo)
      
      cbc%binPath = binPath

      RETURN
      END FUNCTION newCbc
!---------------------------------------------------------------------
      FUNCTION newCbcFl(lst) RESULT(cbc)
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(cbcType) :: cbc

      INTEGER nX
      CHARACTER(LEN=stdL) schm
      TYPE(lstType), POINTER :: lPBC, lPtr
      TYPE(fileType) iniF, binF

      lPBC => lst%get(schm,"Couple to cplBC")
      IF (.NOT.ASSOCIATED(lPBC)) RETURN

      lPtr => lPBC%get(nX,"Number of unknowns",1,ll=0)
      lPtr => lPBC%get(binF,"0D code file path",1)
      lPtr => lPBC%get(iniF,"Unknowns initialization file path",1)

      cbc = newCbc(schm, binF%name(), iniF, nX)

      RETURN
      END FUNCTION newCbcFl
!---------------------------------------------------------------------
!     cbc faces are initialized here
      SUBROUTINE addBcCbc(cbc, name, bGrp)
      CLASS(cbcType), INTENT(INOUT) :: cbc
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: bGrp

      TYPE(cplFaceType), ALLOCATABLE :: cfa(:)

      IF (.NOT.cbc%izd) io%e = "addBcCbc: not initialized"
      cbc%nFa = cbc%nFa + 1
      IF (cbc%nFa .GT. 1) THEN
         IF (ANY(cbc%fa%name .EQ. name)) io%e ="addBcCbc: Repeated name"
         CALL MOVE_ALLOC(cbc%fa,cfa)
         ALLOCATE(cbc%fa(cbc%nFa))
         cbc%fa(1:cbc%nFa-1) = cfa
         DEALLOCATE(cfa)
      ELSE
         ALLOCATE(cbc%fa(cbc%nFa))
      END IF
      cbc%fa(cbc%nFa)%name = name
      cbc%fa(cbc%nFa)%y    = 0D0
      cbc%fa(cbc%nFa)%bGrp = bGrp

      RETURN
      END SUBROUTINE addBcCbc
!---------------------------------------------------------------------
!     true if initialized
      FUNCTION isIzdCbc(cbc) RESULT(izd)
      CLASS(cbcType), INTENT(INOUT) :: cbc
      LOGICAL :: izd

      izd = cbc%izd

      RETURN
      END FUNCTION isIzdCbc
!---------------------------------------------------------------------
!     true if initialized
      FUNCTION isECbc(cbc) RESULT(isE)
      CLASS(cbcType), INTENT(INOUT) :: cbc
      LOGICAL :: isE

      isE = cbc%schm .EQ. cbc_E

      RETURN
      END FUNCTION isECbc
!--------------------------------------------------------------------
!     Solving 0D domain
      SUBROUTINE solveCbc(cbc, U, P)
      CLASS(cbcType), INTENT(INOUT) :: cbc
      TYPE(gVarType), INTENT(IN) :: U, P

      INTEGER iFa, nNo, iM, iCb
      REAL(KIND=8) A
      REAL(KIND=8), ALLOCATABLE :: Uo(:,:), Un(:,:), Po(:), Pn(:),
     2   tmp(:,:)
      TYPE(mshType), POINTER :: lM

      IF (U%dof.NE.nsd .OR. P%dof.NE.1) io%e = 
     2   "solveCbc: Incompatible dof; implementation needed"
      IF (.NOT.ASSOCIATED(U%dmn,P%dmn)) io%e ="solveCbc: U%dmn.NE.P%dmn"

      DO iCb=1, cbc%nFa
         CALL U%dmn%find(cbc%fa(iCb)%name, iM, iFa)
         lM => U%dmn%msh(iM)
         nNo = lM%nNo
         IF (cbc%fa(iCb)%bGrp .EQ. cbc_Neu) THEN
            ALLOCATE(Uo(nsd,nNo), Un(nsd,nNo))
            Uo = U%o%v(:,lM%uP)
            Un = U%n%v(:,lM%uP)
            cbc%fa(iCb)%Qo = lM%integ(iFa,Uo)
            cbc%fa(iCb)%Qn = lM%integ(iFa,Un)
            DEALLOCATE(Uo, Un)
         ELSE IF (cbc%fa(iCb)%bGrp .EQ. cbc_Dir) THEN
            ALLOCATE(Po(nNo), Pn(nNo))
            Po = P%o%s(lM%uP)
            Pn = P%n%s(lM%uP)
            A  = lM%fa(iFa)%area
            cbc%fa(iCb)%Po = lM%integ(iFa,Po)/A
            cbc%fa(iCb)%Pn = lM%integ(iFa,Pn)/A
            DEALLOCATE(Po, Pn)
         END IF
      END DO
      ALLOCATE(tmp(cbc%nFa,4))
      tmp(:,1) = cbc%fa%Qo
      tmp(:,2) = cbc%fa%Qn
      tmp(:,3) = cbc%fa%Po
      tmp(:,4) = cbc%fa%Pn
      tmp = cm%reduce(tmp)
      cbc%fa%Qo = tmp(:,1)
      cbc%fa%Qn = tmp(:,2)
      cbc%fa%Po = tmp(:,3)
      cbc%fa%Pn = tmp(:,4)
      DEALLOCATE(tmp)

      CALL cbc%callLD()
      CALL cbc%calcDer()

      RETURN
      END SUBROUTINE solveCbc
!--------------------------------------------------------------------
!     cplBC derivative is calculated here
      SUBROUTINE calcDerCbc(cbc)
      CLASS(cbcType), INTENT(INOUT) :: cbc

      REAL(KIND=8), PARAMETER :: absTol = 1D-8, relTol = 1D-5
      INTEGER iCb
      REAL(KIND=8) orgQ, orgY, diff
      
      IF (ALL(cbc%fa%bGrp.EQ.cbc_Dir) .OR. cbc%schm.EQ.cbc_E) RETURN
      IF (cbc%schm.EQ.cbc_SI .AND. cbc%called) RETURN

      diff = relTol*NORM2(cbc%fa%Qo)/REAL(cbc%nFa,8)
      IF (diff .LT. absTol) diff = absTol

      DO iCb=1, cbc%nFa
         IF (cbc%fa(iCb)%bGrp .NE. cbc_Neu) CYCLE
         
         orgY = cbc%fa(iCb)%y
         orgQ = cbc%fa(iCb)%Qn
         cbc%fa(iCb)%Qn = cbc%fa(iCb)%Qn + diff

         CALL cbc%callLD()
            
         cbc%fa(iCb)%r  = (cbc%fa(iCb)%y - orgY)/diff
         cbc%fa(iCb)%y  = orgY
         cbc%fa(iCb)%Qn = orgQ
      END DO
      cbc%called = .TRUE.

      RETURN
      END SUBROUTINE calcDerCbc
!--------------------------------------------------------------------
!     Interface to call Lower-D code
      SUBROUTINE callLDCbc(cbc)
      CLASS(cbcType), INTENT(INOUT) :: cbc

      INTEGER iCb, ierr
      CHARACTER(LEN=stdL) ctmp
      TYPE(fileType) f
      REAL(KIND=8), ALLOCATABLE :: y(:)

      ALLOCATE(y(cbc%nFa))
      IF (cm%mas()) THEN
         ctmp = TRIM(appPath)//cbc_cmName
         f = fileType(ctmp,'binary')
         CALL f%open('w')
         WRITE(f%id()) cbc%nFa, cbc%nX, dt
         WRITE(f%id()) cbc%xo
         DO iCb=1, cbc%nFa
            WRITE(f%id()) cbc%fa(iCb)%bGrp, cbc%fa(iCb)%Qo, 
     2         cbc%fa(iCb)%Qn, cbc%fa(iCb)%Po, cbc%fa(iCb)%Pn,
     3         cbc%fa(iCb)%name
         END DO
         CALL f%close()
         CALL SYSTEM(TRIM(cbc%binPath)//" "//TRIM(ctmp))
         CALL f%open('r')
         READ(f%id()) ierr
         IF (ierr .NE. 1) io%e = TRIM(cbc%binPath)//" returned an error"
         READ(f%id()) cbc%xn, y
         CALL f%close()
      END IF
      CALL cm%bcast(cbc%xn)
      CALL cm%bcast(y)
      cbc%fa%y = y

      RETURN
      END SUBROUTINE callLDCbc
!--------------------------------------------------------------------
!     Returns computed values
      SUBROUTINE getCbc(cbc, name, g, r)
      CLASS(cbcType), INTENT(IN) :: cbc
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL(KIND=8), INTENT(OUT) :: g(:), r

      INTEGER iCb

      DO iCb=1, cbc%nFa
         IF (cbc%fa(iCb)%name .EQ. name) THEN
            g = cbc%fa(iCb)%y
            r = cbc%fa(iCb)%r
            EXIT
         END IF
      END DO
      IF (iCb .GT. cbc%nFa) io%e = "getCbc: non-existent fa%name"
      
      RETURN
      END SUBROUTINE getCbc
!---------------------------------------------------------------------
      SUBROUTINE freeCbc(this)
      CLASS(cbcType) :: this

      IF (ALLOCATED(this%fa)) DEALLOCATE(this%fa)
      IF (ALLOCATED(this%xn)) DEALLOCATE(this%xn)
      IF (ALLOCATED(this%xo)) DEALLOCATE(this%xo)
      this%called = .FALSE.
      this%izd    = .FALSE.
      this%nFa    = 0
      this%nX     = 0

      RETURN 
      END SUBROUTINE freeCbc
      END MODULE CBCMOD
