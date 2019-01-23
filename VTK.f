!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Class for generating VTK files. For usage and testing see
!     TEST_VTKMOD at the end of this file. 
!      
!--------------------------------------------------------------------
      MODULE VTKMOD
      USE VARMOD
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: CIntSize = 20

      TYPE, EXTENDS(fileType) :: vtkType
         PRIVATE
!        The pointer to the domain associated with this vtk
         TYPE(dmnType), POINTER :: dmn
      CONTAINS
!        To read/write node-based variables
         PROCEDURE :: nRW => nRWVtk
!        To open a vtk file for either read/write
         PROCEDURE :: open => openVtk
!        Write the partition data
         PROCEDURE :: wPart => wPartVtk
!        To write data of a mesh in raw format
         PROCEDURE, PRIVATE :: wRawRSVtk
         PROCEDURE, PRIVATE :: wRawRVVtk
         GENERIC, PRIVATE :: wRaw => wRawRSVtk, wRawRVVtk
!        To read data of a mesh in raw format
         PROCEDURE, PRIVATE :: rRawRSVtk
         PROCEDURE, PRIVATE :: rRawRVVtk
         GENERIC, PRIVATE :: rRaw => rRawRSVtk, rRawRVVtk
      END TYPE vtkType

      INTERFACE vtkType
         PROCEDURE :: newVtk
      END INTERFACE vtkType

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
      FUNCTION newVtk(name, dmn) RESULT(f)
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      TYPE(vtkType) :: f

      INTEGER i
      CHARACTER(LEN=stdL) stmp

      stmp = ADJUSTL(name)
      i = LEN(TRIM(stmp))
      IF (LOWER(stmp(i-3:i)) .NE. ".vtk") stmp = TRIM(stmp)//'.vtk'

      f%fileType = fileType(TRIM(stmp), 'binary')
      f%dmn => dmn

      RETURN
      END FUNCTION newVtk
!---------------------------------------------------------------------
      SUBROUTINE openVtk(f,rw,name)
      CLASS(vtkType), INTENT(INOUT) :: f
      CHARACTER, INTENT(IN) :: rw
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: name

      INTEGER iM, nEl, nSh, e, cnt, i
      CHARACTER(LEN=stdL) stmp
      INTEGER, ALLOCATABLE :: gIEN(:)
      TYPE(mshType), POINTER :: lM
      
      IF (.NOT.f%dmn%izd) io%e = "openVtk: Un-initialized domain"

      IF (PRESENT(name)) THEN
         stmp = ADJUSTL(name)
         i = LEN(TRIM(stmp))
         IF (stmp(i-3:i) .NE. ".vtk") stmp = TRIM(stmp)//'.vtk'
      ELSE
         stmp = f%name()
      END IF

      CALL f%fileType%open(rw,TRIM(stmp))
      IF (cm%slv()) CALL f%close()
      
!     Reading/writing the header
      CALL f%rw("# vtk DataFile Version 3.0"//eol)
      CALL f%rw("# Prod. by MUPFES"//eol)
      IF (f%isBin()) CALL f%rw("BINARY"//eol)
      IF (.NOT.f%isBin()) CALL f%rw("ASCII"//eol)
      CALL f%rw("DATASET UNSTRUCTURED_GRID"//eol)

!     Reading/writing the Position data
      CALL f%rw("POINTS "//STR(SUM(f%dmn%msh%gnNo),CIntSize)//" float"//
     2   eol)
      DO iM=1, f%dmn%nMsh
         lM => f%dmn%msh(iM)
         IF (f%isRead()) THEN
            CALL f%rRaw(lM%x, lM)
         ELSE
            CALL f%wRaw(lM%x, lM)
         END IF
      END DO

!     No need for slave procs from now on.
      IF (cm%slv()) RETURN

!     Reading/writing the connectivity data
      nEl = SUM(f%dmn%msh%gnEl)
      cnt = SUM((f%dmn%msh%eNoN + 1)*f%dmn%msh%gnEl)
      CALL f%rw("CELLS "//STR(nEl,CIntSize)//" "//STR(cnt,CIntSize)//
     2   eol)
      nSh = -1
      DO iM=1, f%dmn%nMsh
         i  = 1
         lM => f%dmn%msh(iM)
         ALLOCATE(gIEN((lM%eNoN+1)*lM%gnEl))
         IF (.NOT.f%isRead()) THEN
            DO e=1, lM%gnEl
               gIEN(i:i+lM%eNoN) = (/lM%eNoN, lM%gIEN(:,e)+nSh/)
               i = i + lM%eNoN + 1
            END DO
            nSh = nSh + lM%gnNo
         END IF
         CALL f%rw(gIEN)
         DEALLOCATE(gIEN)
      END DO

!     Writing cell types
      CALL f%rw("CELL_TYPES "//STR(nEl,CIntSize)//eol)
      DO iM=1, f%dmn%nMsh
         lM => f%dmn%msh(iM)
         ALLOCATE(gIEN(lM%gnEl))
         IF (.NOT.f%isRead()) gIEN = lM%vtkType
         CALL f%rw(gIEN)
         DEALLOCATE(gIEN)
      END DO
!     Preparing for node-based data write
      CALL f%rw("POINT_DATA "//STR(SUM(f%dmn%msh%gnNo),CIntSize)//eol)

      RETURN
      END SUBROUTINE openVtk
!---------------------------------------------------------------------
      SUBROUTINE wPartVtk(f)
      CLASS(vtkType), INTENT(INOUT) :: f
      
      INTEGER iM, nEl, s, e, i
      INTEGER, ALLOCATABLE :: gP(:)
      TYPE(mshTYpe), POINTER :: lM

      IF (cm%seq() .OR. cm%slv()) RETURN
      IF (.NOT.f%dmn%izd) io%e = "wPartVtk: Un-initialized domain"

      IF (f%isRead() .OR. f%isClose()) 
     2   io%e = "wPartVtk: file must be already open for writing"

      nEl = SUM(f%dmn%msh%gnEl)
      CALL f%rw("CELL_DATA "//STR(nEl,CIntSize)//eol)
      ALLOCATE(gP(nEl))

      CALL f%rw("SCALARS proc-ID int"//eol//"LOOKUP_TABLE default"//eol)
      e = 0
      DO iM=1, f%dmn%nMsh
         lM => f%dmn%msh(iM)
         DO i=0, cm%np() - 1
            s = e + 1
            e = s + lM%eD(i+1) - lM%eD(i) - 1
            gP(s:e) = i
         END DO
      END DO
      CALL f%rw(gP)

      CALL f%rw("SCALARS mesh-ID int"//eol//"LOOKUP_TABLE default"//eol)
      e = 0
      DO iM=1, f%dmn%nMsh
         lM => f%dmn%msh(iM)
         s = e + 1
         e = s + lM%gnEl - 1
         gP(s:e) = iM
      END DO
      CALL f%rw(gP)
      DEALLOCATE(gP)

      RETURN
      END SUBROUTINE wPartVtk
!---------------------------------------------------------------------
      SUBROUTINE nRWVtk(f, var)
      CLASS(vtkType), INTENT(INOUT) :: f
      CLASS(varType), INTENT(INOUT) :: var
 
      INTEGER iM
      REAL(KIND=8), ALLOCATABLE :: rS(:), rV(:,:)
      TYPE(mshTYpe), POINTER :: lM

      IF (.NOT.ASSOCIATED(var%dmn,f%dmn)) io%e = "nRWVtk: var and f"//
     2   " are not associated with the same dmn"

      IF (var%dof .EQ. 1) THEN
         CALL f%rw("SCALARS "//TRIM(var%name)//" float"//eol)
         CALL f%rw("LOOKUP_TABLE default"//eol)
         DO iM=1, f%dmn%nMsh
            lM => f%dmn%msh(iM)
            ALLOCATE(rS(lM%nNo))
            IF (f%isRead()) THEN
               CALL f%rRaw(rS, lM)
               var%s(lM%uP) = rS
            ELSE
               rS = var%s(lM%uP)
               CALL f%wRaw(rS, lM)
            END IF
            DEALLOCATE(rS)
         END DO
      ELSE IF (var%dof .EQ. nsd) THEN
         CALL f%rw("VECTORS "//TRIM(var%name)//" float"//eol)
         DO iM=1, f%dmn%nMsh
            lM => f%dmn%msh(iM)
            ALLOCATE(rV(nsd,lM%nNo))
            IF (f%isRead()) THEN
               CALL f%rRaw(rV, lM)
               var%v(:,lM%uP) = rV
            ELSE
               rV = var%v(:,lM%uP)
               CALL f%wRaw(rV, lM)
            END IF
            DEALLOCATE(rV)
         END DO
      ELSE
         io%e = "nRWVtk: unknown var%dof"
      END IF

      RETURN
      END SUBROUTINE nRWVtk

!#####################################################################

      SUBROUTINE wRawRSVtk(f, r, lM)
      CLASS(vtkType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(IN) :: r(:)
      TYPE(mshType), INTENT(IN) :: lM
 
      LOGICAL flag
      INTEGER Ac
      REAL(KIND=8), ALLOCATABLE :: gS(:)
      REAL, ALLOCATABLE :: gSR(:)

      ALLOCATE(gS(lM%gnNo))
      CALL cm%global(r, gS, lM%gN)
      IF (cm%slv()) RETURN

      ALLOCATE(gSR(SIZE(gS)))
      gSR = REAL(gS)

      flag = .TRUE.
      DO Ac=1, SIZE(gSR)
         IF (gSR(Ac).NE.gSR(Ac) .OR. gSR(Ac).GT.HUGE(gSR)) THEN
            IF (flag) THEN
               io%w = "Exceptional number detected, replaced by zero"
               flag = .FALSE.
            END IF
            gSR(Ac) = 0D0
         END IF
      END DO
      
      CALL f%rw(gSR)
      DEALLOCATE(gS,gSR)

      RETURN
      END SUBROUTINE wRawRSVtk
!---------------------------------------------------------------------
      SUBROUTINE rRawRSVtk(f, r, lM)
      CLASS(vtkType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(INOUT) :: r(:)
      TYPE(mshType), INTENT(IN) :: lM
 
      REAL(KIND=8), ALLOCATABLE :: gS(:)
      REAL, ALLOCATABLE :: gSR(:)

      IF (cm%mas()) THEN
         ALLOCATE(gSR(lM%gnNo))
         CALL f%rw(gSR)
         ALLOCATE(gS(lM%gnNo))
         gS = REAL(gSR,8)
         DEALLOCATE(gSR)
      ELSE
         ALLOCATE(gS(0))
      END IF

      CALL cm%local(r, gS, lM%gN)
      DEALLOCATE(gS)

      RETURN
      END SUBROUTINE rRawRSVtk
!---------------------------------------------------------------------
      SUBROUTINE wRawRVVtk(f, r, lM)
      CLASS(vtkType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(IN) :: r(:,:)
      TYPE(mshType), INTENT(IN) :: lM
 
      LOGICAL flag
      INTEGER Ac, i
      REAL(KIND=8), ALLOCATABLE :: gV(:,:)
      REAL, ALLOCATABLE :: gSR(:)

      ALLOCATE(gV(SIZE(r,1),lM%gnNo))
      CALL cm%global(r, gV, lM%gN)
      IF (cm%slv()) RETURN

      ALLOCATE(gSR(lM%gnNo*3))
      gSR = 0D0
      DO Ac=1, lM%gnNo
         DO i=1, nsd
            gSR((Ac-1)*3 + i) = REAL(gV(i,Ac))
         END DO
      END DO

      flag = .TRUE.
      DO Ac=1, SIZE(gSR)
         IF (gSR(Ac).NE.gSR(Ac) .OR. gSR(Ac).GT.HUGE(gSR)) THEN
            IF (flag) THEN
               io%w = "Exceptional number detected, replaced by zero"
               flag = .FALSE.
            END IF
            gSR(Ac) = 0D0
         END IF
      END DO
      
      CALL f%rw(gSR)
      DEALLOCATE(gV,gSR)

      RETURN
      END SUBROUTINE wRawRVVtk
!---------------------------------------------------------------------
      SUBROUTINE rRawRVVtk(f, r, lM)
      CLASS(vtkType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(INOUT) :: r(:,:)
      TYPE(mshType), INTENT(IN) :: lM
 
      INTEGER i, Ac
      REAL(KIND=8), ALLOCATABLE :: gV(:,:)
      REAL, ALLOCATABLE :: gSR(:)

      IF (cm%mas()) THEN
         ALLOCATE(gSR(3*lM%gnNo))
         CALL f%rw(gSR)
         ALLOCATE(gV(nsd,lM%gnNo))
         DO Ac=1, lM%gnNo
            DO i=1, nsd
               gV(i,Ac) = REAL(gSR(3*(Ac-1) + i),8)
            END DO
         END DO
         DEALLOCATE(gSR)
      ELSE
         ALLOCATE(gV(0,0))
      END IF

      CALL cm%local(r, gV, lM%gN)
      DEALLOCATE(gV)

      RETURN
      END SUBROUTINE rRawRVVtk
!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_VTKMOD(prs, vel)
      TYPE(varType), INTENT(INOUT) :: prs, vel

      TYPE(vtkType) :: vtk

      prs%s = 1D0
      vel%v = 2D0
      vtk = vtkType('hello',prs%dmn)
      io%o = "newVtk: "//CLR("(PASSED)",3)
      CALL vtk%open('w')
      io%o = "vtk%write-open: "//CLR("(PASSED)",3)
      CALL vtk%nRW(prs)
      CALL vtk%nRW(vel)
      io%o = "vtk%write-nRW: "//CLR("(PASSED)",3)
      CALL vtk%wPart()
      io%o = "vtk%wPart: "//CLR("(PASSED)",3)
      CALL vtk%close()
      io%o = "vtk%close: "//CLR("(PASSED)",3)
      CALL vtk%open('r')
      io%o = "vtk%read-open: "//CLR("(PASSED)",3)
      prs%s = 0D0
      vel%v = 0D0
      CALL vtk%nRW(prs)
      CALL vtk%nRW(vel)
      IF (.NOT.ISZERO(prs%sP,prs%sP) .OR.
     2    .NOT.ISZERO(vel%sP,vel%sP)) io%e = "Issue with vtk%nRW"
      io%o = "vtk%read-nRW: "//CLR("(PASSED)",3)
      CALL vtk%close()
      io%o = "Test completed. Also inspect hello.vtk."
      
      RETURN
      END SUBROUTINE TEST_VTKMOD
      END MODULE VTKMOD
