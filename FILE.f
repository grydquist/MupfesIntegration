!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     File class for writing and reading. Handles any permutation of 
!     ASCII/Binary and read/write. For usage and testing the
!     class see TEST_FILEMOD at the end of this file. 
!      
!--------------------------------------------------------------------
      MODULE FILEMOD
      USE UTILMOD
      IMPLICIT NONE

      TYPE fileType
         PRIVATE
!        Whether the file is declared yet
         LOGICAL :: izd = .FALSE.
!        Whether the file is opened yet
         LOGICAL :: isO = .FALSE.
!        Whether the file is binary
         LOGICAL :: isB = .FALSE.
!        Whether writing
         LOGICAL :: isW = .FALSE.
!        File ID
         INTEGER :: fid = -1
!        Positin of last written data (for bindary)
         INTEGER :: pos = 1
!        File name
         CHARACTER(LEN=stdL) :: fName = DEFAULT_NAME
      CONTAINS 
!        opens the file
         PROCEDURE :: open => openFile
!        closes the file
         PROCEDURE :: close => closeFile
!        Reads or write depending on how the file is opened
         PROCEDURE :: rwISFile
         PROCEDURE :: rwRSFile
         PROCEDURE :: rwIVFile
         PROCEDURE :: rwRVFile
         PROCEDURE :: rwDVFile
         PROCEDURE :: rwDMFile
         PROCEDURE :: rwSSFile
         GENERIC :: rw => rwISFile, rwRSFile, rwIVFile, rwRVFile,
     2      rwDVFile, rwDMFile, rwSSFile
!        Returns file ID
         PROCEDURE :: id => idFile
!        Returns true if file is close
         PROCEDURE :: isClose => isCloseFile
!        Returns true if file is binary
         PROCEDURE :: isBin => isBinFile
!        Returns true if is opened for read
         PROCEDURE :: isRead => isReadFile
!        Returns true if the file exists
         PROCEDURE :: exist => existFile
!        Returns file name
         PROCEDURE :: name => nameFile
!        Set the current position (for binary files)
         PROCEDURE :: setPos => setPosFile
      END TYPE fileType

      INTERFACE fileType
         PROCEDURE :: newFile
      END INTERFACE fileType

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
      FUNCTION newFile(name, format) RESULT(f)
      CHARACTER(LEN=*), INTENT(IN) :: name
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: format
      TYPE(fileType) :: f

      INTEGER i
      CHARACTER(LEN=stdL) stmp

      f%izd   = .TRUE.
      f%fName = ADJUSTL(name)
      IF (PRESENT(format)) THEN
         SELECT CASE (LOWER(format))
         CASE ('binary')
            f%isB = .TRUE.
         CASE ('ascii')
            f%isB = .FALSE.
         CASE DEFAULT
            ERROR STOP "newFile: Unknown format"
         END SELECT
      ELSE
         stmp = f%name()
         i = LEN(TRIM(stmp))
         IF (i .GT. 3) THEN
            f%isB = LOWER(stmp(i-3:i)).EQ.'.bin'
         END IF
      END IF

      RETURN
      END FUNCTION newFile
!####################################################################
      SUBROUTINE openFile(f,rw,name) 
      CLASS(fileType), INTENT(INOUT) :: f
      CHARACTER, INTENT(IN) :: rw
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: name

      LOGICAL flag
      INTEGER fid

!     If file is not declared, we return without doing anything
      IF (.NOT.f%izd) RETURN
      
      IF (PRESENT(name)) f%fName = ADJUSTL(name)

      IF (LOWER(rw) .EQ. 'r') THEN
         f%isW = .FALSE.
      ELSE IF (LOWER(rw).EQ.'w' .OR. LOWER(rw).EQ.'a') THEN
         f%isW = .TRUE.
      ELSE
         ERROR STOP "openFile: acceptable rw are 'r' and 'w'"
      END IF
      IF (.NOT.f%isClose()) THEN
         PRINT *, "openFile: file <"//TRIM(f%fName)//"> is already open"
         RETURN
      END IF
      f%isO = .TRUE.

      IF (.NOT.f%exist() .AND. .NOT.f%isW) THEN
         PRINT *, "openFile: file <"//TRIM(f%fName)//
     2      "> does not exists or can not be opened"
         ERROR STOP
      END IF

      DO fid=11, 1024
         INQUIRE(UNIT=fid, OPENED=flag)
         IF (.NOT.flag) EXIT
      END DO
      f%fid = fid
      IF (f%isBin()) THEN
         IF (LOWER(rw) .EQ. 'w') THEN
            OPEN(UNIT=f%fid, FILE=TRIM(f%fName))
            WRITE(f%fid,*) ""
            CLOSE(f%fid)
         END IF
         IF (LOWER(rw) .EQ. 'a') THEN
            OPEN(UNIT=f%fid, FILE=TRIM(f%fName), ACCESS="STREAM", 
     2         CONVERT='BIG_ENDIAN', POSITION='APPEND')
            INQUIRE(UNIT=f%fid,POS=f%pos)
         ELSE
            OPEN(UNIT=f%fid, FILE=TRIM(f%fName), ACCESS="STREAM", 
     2         CONVERT='BIG_ENDIAN')
         END IF
      ELSE
         IF (LOWER(rw) .EQ. 'a') THEN
            OPEN(UNIT=f%fid, FILE=TRIM(f%fName), POSITION='APPEND')
         ELSE
            OPEN(UNIT=f%fid, FILE=TRIM(f%fName))
         END IF
      END IF

      RETURN
      END SUBROUTINE openFile
!---------------------------------------------------------------------
      SUBROUTINE closeFile(f)
      CLASS(fileType), INTENT(INOUT) :: f

      IF (f%isClose()) RETURN

      CALL f%setPos(1)
      CLOSE(f%fid)
      f%fid = -1
      f%isO = .FALSE.

      RETURN
      END SUBROUTINE closeFile
!---------------------------------------------------------------------
      PURE FUNCTION idFile(f) RESULT(fid)
      CLASS(fileType), INTENT(IN) :: f
      INTEGER fid

      fid = f%fid

      RETURN
      END FUNCTION idFile
!---------------------------------------------------------------------
      PURE FUNCTION isBinFile(f) RESULT(isB)
      CLASS(fileType), INTENT(IN) :: f
      LOGICAL isB

      isB = f%isB

      RETURN
      END FUNCTION isBinFile
!---------------------------------------------------------------------
      PURE FUNCTION isCloseFile(f) RESULT(isC)
      CLASS(fileType), INTENT(IN) :: f
      LOGICAL isC

      isC = .NOT.f%isO

      RETURN
      END FUNCTION isCloseFile
!---------------------------------------------------------------------
      PURE FUNCTION isReadFile(f) RESULT(isR)
      CLASS(fileType), INTENT(IN) :: f
      LOGICAL isR

      isR = .NOT.f%isW

      RETURN
      END FUNCTION isReadFile
!---------------------------------------------------------------------
      FUNCTION existFile(f) RESULT(flag)
      CLASS(fileType), INTENT(IN) :: f
      LOGICAL flag

      INQUIRE(FILE=TRIM(f%fName), EXIST=flag)

      RETURN
      END FUNCTION existFile
!---------------------------------------------------------------------
      PURE FUNCTION nameFile(f) RESULT(fName)
      CLASS(fileType), INTENT(IN) :: f
      CHARACTER(LEN=stdL) fName

      fName = f%fName

      RETURN
      END FUNCTION nameFile
!---------------------------------------------------------------------
      SUBROUTINE setPosFile(f, pos)
      CLASS(fileType), INTENT(INOUT) :: f
      INTEGER, INTENT(IN) :: pos

      f%pos = pos

      RETURN
      END SUBROUTINE setPosFile
!---------------------------------------------------------------------
      SUBROUTINE rwISFile(f,r)
      CLASS(fileType), INTENT(INOUT) :: f
      INTEGER, INTENT(INOUT) :: r

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) r
         ELSE
            READ(f%fid,POS=f%pos) r
         END IF
         f%pos = f%pos + SIZEOF(r)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,*) r
         ELSE
            READ(f%fid,*) r
         END IF
      END IF

      RETURN
      END SUBROUTINE rwISFile
!---------------------------------------------------------------------
      SUBROUTINE rwRSFile(f,r)
      CLASS(fileType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(INOUT) :: r

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) r
         ELSE
            READ(f%fid,POS=f%pos) r
         END IF
         f%pos = f%pos + SIZEOF(r)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,*) r
         ELSE
            READ(f%fid,*) r
         END IF
      END IF

      RETURN
      END SUBROUTINE rwRSFile
!---------------------------------------------------------------------
      SUBROUTINE rwIVFile(f,r)
      CLASS(fileType), INTENT(INOUT) :: f
      INTEGER, INTENT(INOUT) :: r(:)

      INTEGER rSpl

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) r
         ELSE
            READ(f%fid,POS=f%pos) r
         END IF
         f%pos = f%pos + SIZEOF(rSpl)*SIZE(r)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,*) r
         ELSE
            READ(f%fid,*) r
         END IF
      END IF

      RETURN
      END SUBROUTINE rwIVFile
!---------------------------------------------------------------------
      SUBROUTINE rwRVFile(f,r)
      CLASS(fileType), INTENT(INOUT) :: f
      REAL, INTENT(INOUT) :: r(:)

      REAL rSpl

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) r
         ELSE
            READ(f%fid,POS=f%pos) r
         END IF
         f%pos = f%pos + SIZEOF(rSpl)*SIZE(r)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,*) r
         ELSE
            READ(f%fid,*) r
         END IF
      END IF

      RETURN
      END SUBROUTINE rwRVFile
!---------------------------------------------------------------------
      SUBROUTINE rwDVFile(f,r)
      CLASS(fileType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(INOUT) :: r(:)

      REAL(KIND=8) rSpl

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) r
         ELSE
            READ(f%fid,POS=f%pos) r
         END IF
         f%pos = f%pos + SIZEOF(rSpl)*SIZE(r)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,*) r
         ELSE
            READ(f%fid,*) r
         END IF
      END IF

      RETURN
      END SUBROUTINE rwDVFile
!---------------------------------------------------------------------
      SUBROUTINE rwDMFile(f,r)
      CLASS(fileType), INTENT(INOUT) :: f
      REAL(KIND=8), INTENT(INOUT) :: r(:,:)

      REAL(KIND=8) rSpl

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) r
         ELSE
            READ(f%fid,POS=f%pos) r
         END IF
         f%pos = f%pos + SIZEOF(rSpl)*SIZE(r)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,*) r
         ELSE
            READ(f%fid,*) r
         END IF
      END IF

      RETURN
      END SUBROUTINE rwDMFile
!---------------------------------------------------------------------
      SUBROUTINE rwSSFile(f,s,r)
      CLASS(fileType), INTENT(INOUT) :: f
      CHARACTER(LEN=*), INTENT(IN) :: s
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: r

      CHARACTER(LEN=LEN(s)) fr

      IF (f%isClose()) RETURN
      IF (f%isBin()) THEN
         IF (f%isW) THEN
            WRITE(f%fid,POS=f%pos) s
         ELSE
            READ(f%fid,POS=f%pos) fr
         END IF
         f%pos = f%pos + LEN(s)*SIZEOF(eol)
      ELSE
         IF (f%isW) THEN
            WRITE(f%fid,'(A)') s
         ELSE
            READ(f%fid,'(A)') fr
         END IF
      END IF
      IF (.NOT.f%isW) THEN
         IF (PRESENT(r)) THEN
            r = fr
         ELSE
            IF (fr .NE. s) THEN
               print *, len(s), len(fr)
               PRINT *, "<"//TRIM(fr)//"> .NE. <"//TRIM(s)//">"
               ERROR STOP "rwSSFile-R: Non-matching strings"
            END IF
         END IF
      END IF

      RETURN
      END SUBROUTINE rwSSFile
!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_FILEMOD()

      TYPE(fileType) fA, fB
      REAL(KIND=8) :: r(16)

      fA = fileType('Hello-A.txt') ! ASCII
      fB = fileType('Hello-B.txt','binary') ! Binary
      IF (fA%name().NE.'Hello-A.txt') ERROR STOP "Issue with newFile"
      PRINT *, "newFile: "//CLR("(PASSED)",3)
      CALL fA%open('w')
      CALL fB%open('w')
      PRINT *, "file%open: "//CLR("(PASSED)",3)
      r = 1D0
      CALL fA%rw(r)
      CALL fB%rw(r)
      PRINT *, "file%write: "//CLR("(PASSED)",3)
      CALL fA%close()
      CALL fB%close()
      PRINT *, "file%close: "//CLR("(PASSED)",3)
      CALL fA%open('r')
      CALL fB%open('r')
      r = 0D0
      CALL fA%rw(r)
      IF (ANY(r .NE. 1D0)) ERROR STOP "Issue with %read-A"
      r = 0D0
      CALL fB%rw(r)
      IF (ANY(r .NE. 1D0)) ERROR STOP "Issue with %read-B"
      CALL fA%close()
      CALL fB%close()
      PRINT *, "file%read: "//CLR("(PASSED)",3)
!     Cleaning up after myself
      OPEN(1,FILE='Hello-A.txt')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='Hello-B.txt')
      CLOSE(1,STATUS='DELETE')

      RETURN
      END SUBROUTINE TEST_FILEMOD
      END MODULE FILEMOD
