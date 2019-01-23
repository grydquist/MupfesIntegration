!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     This module is used for I/O. For usage see TEST_LSTMOD at the end
!     of this file. 
!      
!--------------------------------------------------------------------

      MODULE LSTMOD
      USE IOMOD
      USE FILEMOD
      USE CMMOD
      IMPLICIT NONE

      TYPE lstType
         PRIVATE
!        Whether this line has been used sofar
         LOGICAL :: used = .FALSE.
!        Length of sub lst
         INTEGER :: l = 0
!        Line number associated with this lst
         INTEGER line
!        Command 
         CHARACTER(LEN=stdL) :: kwd = 'NONE'
!        Value associated with the command
         CHARACTER(LEN=stdL) :: val = 'NONE'
!        Sublst, under current lst
         TYPE(lstType), POINTER :: sub(:)
      CONTAINS 
         PROCEDURE :: ping => pingLst
         PROCEDURE, PRIVATE :: lsrch => lsrchLst
         PROCEDURE :: srch => srchLst
         PROCEDURE :: check => checkLst
         PROCEDURE :: GFLL
         PROCEDURE :: GFLI
         PROCEDURE :: GFLR
         PROCEDURE :: GFLS
         PROCEDURE :: GFLF
         PROCEDURE :: GFLV
         GENERIC :: get => GFLL, GFLI, GFLR, GFLS, GFLF, GFLV
      END TYPE lstType
      
      INTERFACE lstType
         PROCEDURE :: newLst
      END INTERFACE lstType

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
      FUNCTION newLst(fileName) RESULT(lst)
      CHARACTER(LEN=*), INTENT(IN) :: fileName
      TYPE(lstType) lst
      
      INTEGER, PARAMETER :: maxL = 100

      LOGICAL flag
      INTEGER s, e, l, i, j, A, fNl, lInd(stdL), lvl, fid
      CHARACTER(LEN=stdL) ctmp, sTmp
      
      INTEGER, ALLOCATABLE :: fInd(:), fLineN(:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: fCon(:)

      lst%kwd  = fileName
      lst%used = .TRUE.
 
      IF (cm%mas()) THEN
!     These are the default values used for the first few lines below
         fid = 1
         INQUIRE(FILE=TRIM(fileName), EXIST=flag)
         IF (.NOT.flag) THEN
            io%e = "File <"//TRIM(fileName)//"> does not exist"//
     2         " or can not be opened"
         END IF
         OPEN(fid, FILE=fileName, STATUS="OLD")

!     Reading the entire ini file      
         fNl = 0
         DO
            READ(fid,"(A)",IOSTAT=i) sTmp
            sTmp = ADJUSTC(sTmp)
            IF (i .GT. 0) THEN
               io%e = "While reading: "//fileName
            ELSE IF (i .LT. 0) THEN
               EXIT
            END IF
            IF (sTmp .EQ. "End") EXIT
            fNl = fNl + 1
         END DO

!     This for the case a curly brackets is at the middle of a line
         e = 100 
         ALLOCATE(fCon(fNl+e), fInd(fNl+e), fLineN(fNl+e))
         fCon = ""
         REWIND(fid)
         A = 0
         lvl = 1
         DO i=1, fNl
            READ(fid,"(A)") sTmp
            sTmp = ADJUSTC(sTmp)
            l = LEN(TRIM(sTmp))
!     Empty and commented lines are ignored
            IF (l.EQ.0 .OR. sTmp(1:1).EQ."#") CYCLE
         
            e = 0
            ctmp = ""
!     Going one level up or down if "{" or "}" are found
            DO j=1, l
               IF (sTmp(j:j) .EQ. "{") THEN
                  lvl = lvl + 1
               ELSE IF (sTmp(j:j) .EQ. "}") THEN
                  lvl = lvl - 1
               ELSE
                  e = e + 1
                  lInd(e) = lvl
                  ctmp(e:e) = sTmp(j:j)
               END IF
            END DO
            l = e
            sTmp = ctmp
            IF (l .EQ. 0) CYCLE
!     Cutting a line in pieces if the level is changing in a particular
!     character in the line
            e = lInd(1)
            s = 1
            DO j=1, l
               IF (lInd(j) .NE. e) THEN
                  A = A + 1
                  fCon(A) = ADJUSTC(sTmp(s:j-1))
                  fInd(A) = e
                  fLineN(A) = i
                  s = j
                  e = lInd(j)
               END IF
            END DO
            A = A + 1
            fCon(A) = ADJUSTC(sTmp(s:l))
            fInd(A) = e
            fLineN(A) = i
         END DO
         IF (ANY(fInd(1:A) .LT. 1)) THEN
            io%e = "Too many } in "//fileName
         END IF
         IF (lvl .NE. 1) io%e = "Too many { in "//fileName
         CLOSE(fid)
         fNl = A
         i   = SIZE(fInd)
      END IF
      CALL cm%bcast(fNl)
      CALL cm%bcast(i)
      IF (cm%slv()) ALLOCATE(fCon(i), fInd(i), fLineN(i))
      CALL cm%bcast(fCon)
      CALL cm%bcast(fInd)
      CALL cm%bcast(fLineN)

!     transforming the entire thing into a lst
      CALL SETSUBLST(1, fNl, lst, fCon, fInd, fLineN)
      DEALLOCATE(fCon, fInd, fLineN)
      
      RETURN
      END FUNCTION newLst
!--------------------------------------------------------------------
!     This subroutine sets the length of a lst, allocates and sets 
!     value to lst%sub based on the lines that are in the same 
!     level as first line. To call this, you need to provide the lst 
!     structure, first line number (s) and last line number (e) 
      RECURSIVE SUBROUTINE SETSUBLST(s, e, lst, fCon, fInd, fLineN)
      INTEGER, INTENT(IN) :: s, e
      CHARACTER(LEN=stdL), INTENT(IN) :: fCon(:)
      INTEGER, INTENT(IN) :: fInd(:), fLineN(:)
      TYPE(lstType), INTENT(INOUT) :: lst

      LOGICAL flag
      INTEGER lvl, i, j, k, l
      CHARACTER(LEN=stdL) sTmp

!     Counting the numebr of lines to be allocated
      lvl = fInd(s)
      DO i=s, e
         IF (fInd(i) .EQ. lvl) lst%l = lst%l + 1
      END DO

      ALLOCATE(lst%sub(lst%l))
      j = 0
!     We assume this lst does not have a sublst.      
      flag = .FALSE.
      DO i=s, e
         IF (fInd(i) .EQ. lvl) THEN
            IF (flag) THEN
!     This is the end of sublst               
               flag = .FALSE.
               CALL SETSUBLST(k, i-1, lst%sub(j), fCon, fInd, fLineN)
            END IF
            j = j + 1
            sTmp = fCon(i)
            DO l=1, LEN(TRIM(sTmp))
               IF (sTmp(l:l) .EQ. ":") EXIT
            END DO
            lst%sub(j)%kwd  = LOWER(sTmp(:l-1))
            lst%sub(j)%val  = ADJUSTC(sTmp(l+1:))
            lst%sub(j)%line = fLineN(i)
         ELSE
!     There is a sublst!            
            IF (.NOT. flag) k = i
            flag = .TRUE.
         END IF
      END DO
!     This is for the case that index of last line, fInd(e), is
!     different from lvl  
      IF (flag) CALL SETSUBLST(k, i-1, lst%sub(j), fCon, fInd, fLineN)

      RETURN
      END SUBROUTINE SETSUBLST
!####################################################################
!     This function search through lst for "cmnd" and returns the line
!     number that is found. The index of the searched line can be 
!     inputted by iInd. If iInd is provided, a line is not found, error
!     is thrown
      FUNCTION lsrchLst(lst, cmnd, iInd) RESULT(lsrch)
      CLASS(lstType), INTENT(INOUT) :: lst
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: iInd
      INTEGER lsrch
      
      INTEGER i, n, ind
      
      ind = 1
      IF (PRESENT(iInd)) ind = iInd

!     If the line is not found, -1 is returned
      lsrch = -1
      n = 0
      DO i=1, lst%l
         IF (lst%sub(i)%kwd .EQ. LOWER(cmnd)) THEN
            n = n + 1
            IF (n .EQ. ind) THEN
               lst%sub(i)%used = .TRUE.
               lsrch = i
               RETURN
            END IF
         END IF
      END DO
      IF (PRESENT(iInd)) THEN
         io%e = TRIM(lst%ping(cmnd))//" Command not found"
      END IF

      RETURN
      END FUNCTION lsrchLst
!--------------------------------------------------------------------
      FUNCTION srchLst(lst, cmnd, ll, ul) RESULT(srch)
      CLASS(lstType), INTENT(IN) :: lst
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ll, ul
      INTEGER srch
      
      INTEGER i

      srch = 0
      DO i=1, lst%l
         IF (lst%sub(i)%kwd .EQ. LOWER(cmnd)) srch = srch + 1
      END DO

      IF (PRESENT(ll)) THEN
         IF (srch .LT. ll) io%e = TRIM(lst%ping(cmnd))//
     2      " COMMAND must be repeated more/equal than "//ll//" time/s"
      END IF
      IF (PRESENT(ul)) THEN
         IF (srch .GT. ul) io%e = TRIM(lst%ping(cmnd))//
     2      " COMMAND must be repeated less/equal than "//ul//" time/s"
      END IF
      
      RETURN
      END FUNCTION srchLst

!####################################################################
      FUNCTION pingLst(lst, cmnd, subL) RESULT(ping)
      CLASS(lstType), INTENT(IN) :: lst
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      TYPE(lstType), OPTIONAL :: subL
      CHARACTER(LEN=stdL) ping

      IF (PRESENT(subL)) THEN
         ping = "At LINE "//subL%line//", searched COMMAND <"//
     2      TRIM(cmnd)//"> under <"//TRIM(lst%kwd)//": "//
     3      TRIM(lst%val)//"> ::"
      ELSE
         ping = "Near LINE "//lst%line//" at <"//TRIM(lst%kwd)//
     2      ": "//TRIM(lst%val)//">, searched COMMAND <"//
     3      TRIM(cmnd)//"> ::"
      END IF
 
      RETURN
      END FUNCTION pingLst

!####################################################################
!     This routine checks all the lines of the lst to make sure they
!     have been used
      RECURSIVE SUBROUTINE checkLst(lst)
      CLASS(lstType), INTENT(IN) :: lst

      INTEGER i

      IF (.NOT.lst%used) io%e = "Near line "//lst%line//
     2   " keyword <"//TRIM(lst%kwd)//"> is not recognized"

      DO i=1, lst%l
         CALL checkLst(lst%sub(i))
      END DO

      RETURN
      END SUBROUTINE checkLst

!####################################################################
!     This is for reading logical
      FUNCTION GFLL(lst, lVal, cmnd, ind)
      CLASS(lstType), INTENT(INOUT) :: lst
      LOGICAL, INTENT(OUT) :: lVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(lstType), POINTER :: GFLL

      INTEGER i
      CHARACTER(LEN=stdL) c

      IF (PRESENT(ind)) THEN
         i = lst%LSRCH(cmnd, ind)
         GFLL => lst%GFLS(c, cmnd, ind)
      ELSE
         i = lst%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLL => NULL()
            RETURN
         END IF
         GFLL => lst%GFLS(c, cmnd, 1)
      END IF
      
      IF (c.EQ."1" .OR. c.EQ.'true' .OR. c.EQ.'t') THEN
         lVal = .TRUE.
         io%d = TRIM(lst%ping(cmnd,GFLL))//" Read TRUE"
      ELSE IF (c.EQ."0" .OR. c.EQ.'false' .OR. c.EQ.'f') THEN
         lVal = .FALSE.
         io%d = TRIM(lst%ping(cmnd,GFLL))//" Read FALSE"
      ELSE
         io%e = TRIM(lst%ping(cmnd,GFLL))//" Reading error"
      END IF
      
      RETURN
      END FUNCTION GFLL
!--------------------------------------------------------------------
!     This function return the value of lstPtr at a particular line 
!     (ind) if ll/lb or ul/ub are specified, the value is checked 
!     to be in that specified range. This is for reading integers.
      FUNCTION GFLI(lst, iVal, cmnd, ind, ll, ul)
      CLASS(lstType), INTENT(INOUT) :: lst
      INTEGER, INTENT(OUT) :: iVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind, ll, ul
      TYPE(lstType), POINTER :: GFLI

      INTEGER i, ioS
 
      IF (PRESENT(ind)) THEN
         i = lst%LSRCH(cmnd, ind)
      ELSE
         i = lst%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLI => NULL()
            RETURN
         END IF
      END IF
      GFLI => lst%sub(i)

      READ(GFLI%val, *, IOSTAT=ioS) iVal
      IF (ioS .NE. 0) THEN
         io%e = TRIM(lst%ping(cmnd,GFLI))//" Reading error"
      END IF
      IF (PRESENT(ll)) THEN
         IF (iVal .LT. ll) io%e = 
     2      TRIM(lst%ping(cmnd,GFLI))//" Lower limit is "//ll
      END IF
      IF (PRESENT(ul)) THEN
         IF (iVal .GT. ul) io%e =
     2      TRIM(lst%ping(cmnd,GFLI))//" Upper limit is "//ul
      END IF
      io%d = TRIM(lst%ping(cmnd,GFLI))//" Read integer value "//
     2   iVal

      RETURN
      END FUNCTION GFLI
!--------------------------------------------------------------------
!     This is for reading real numbers
      FUNCTION GFLR(lst, rVal, cmnd, ind, ll, ul, lb, ub)
      CLASS(lstType), INTENT(INOUT) :: lst
      REAL(KIND=8), INTENT(OUT) :: rVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      REAL(KIND=8), INTENT(IN), OPTIONAL :: ll, ul, lb, ub
      TYPE(lstType), POINTER :: GFLR

      INTEGER i, ioS
 
      IF (PRESENT(ind)) THEN
         i = lst%LSRCH(cmnd, ind)
      ELSE
         i = lst%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLR => NULL()
            RETURN
         END IF
      END IF
      GFLR => lst%sub(i)

      READ(GFLR%val, *, IOSTAT=ioS) rVal
      IF (ioS .NE. 0) THEN
         io%e = TRIM(lst%ping(cmnd,GFLR))//" Reading error"
      END IF
      IF (PRESENT(ll)) THEN
         IF (rVal .LT. ll) io%e = 
     2      TRIM(lst%ping(cmnd,GFLR))//" Lower limit is "//ll
      END IF
      IF (PRESENT(ul)) THEN
         IF (rVal .GT. ul) io%e =
     2      TRIM(lst%ping(cmnd,GFLR))//" Upper limit is "//ul
      END IF
      IF (PRESENT(lb)) THEN
         IF (rVal .LE. lb) io%e = 
     2      TRIM(lst%ping(cmnd,GFLR))//" Lower bound is "//lb
      END IF
      IF (PRESENT(ub)) THEN
         IF (rVal .GE. ub) io%e =
     2      TRIM(lst%ping(cmnd,GFLR))//" Upper bound is "//ub
      END IF
      io%d = TRIM(lst%ping(cmnd,GFLR))//" Read real value "//
     2   rVal
         
      RETURN
      END FUNCTION GFLR
!--------------------------------------------------------------------
!     This is for reading strings      
      FUNCTION GFLS(lst, sVal, cmnd, ind)
      CLASS(lstType), INTENT(INOUT) :: lst
      CHARACTER(LEN=stdL), INTENT(OUT) :: sVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(lstType), POINTER :: GFLS

      INTEGER i, ioS

      IF (PRESENT(ind)) THEN
         i = lst%LSRCH(cmnd, ind)
      ELSE
         i = lst%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLS => NULL()
            RETURN
         END IF
      END IF
      GFLS => lst%sub(i)
      
      READ(GFLS%val,*,IOSTAT=ioS) sVal
      sVal = LOWER(sVal)
      IF (ioS .NE. 0) THEN
         io%e = TRIM(lst%ping(cmnd,GFLS))//" Reading error"
      END IF
      io%d = TRIM(lst%ping(cmnd,GFLS))//" Read char value <"
     2   //TRIM(sVal)//">"

      RETURN
      END FUNCTION GFLS
!--------------------------------------------------------------------
!     This is for reading a vector. The results are returned into v
      FUNCTION GFLV(lst, vVal, cmnd, ind)
      CLASS(lstType), INTENT(INOUT) :: lst
      REAL(KIND=8), INTENT(OUT) :: vVal(:)
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(lstType), POINTER :: GFLV

      INTEGER i, ioS, n
 
      n = SIZE(vVal)
      IF (PRESENT(ind)) THEN
         i = lst%LSRCH(cmnd, ind)
      ELSE
         i = lst%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLV => NULL()
            RETURN
         END IF
      END IF
      GFLV => lst%sub(i)

      READ(GFLV%val, *, IOSTAT=ioS) vVal
      IF (ioS .NE. 0) THEN
         io%e = TRIM(lst%ping(cmnd,GFLV))//" Reading error"
      END IF
      DO i=1, n
         io%d = TRIM(lst%ping(cmnd,GFLV))//" Read vector value "//
     2      vVal(i)
      END DO

      RETURN
      END FUNCTION GFLV
!--------------------------------------------------------------------
!     This function opens a file and returns the handle to the file
      FUNCTION GFLF(lst, fVal, cmnd, ind)
      CLASS(lstType), INTENT(INOUT) :: lst
      TYPE(fileType), INTENT(INOUT) :: fVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(lstType), POINTER :: GFLF

      INTEGER i
      CHARACTER(LEN=stdL) fName

      IF (PRESENT(ind)) THEN
         i = lst%LSRCH(cmnd, ind)
      ELSE
         i = lst%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLF => NULL()
            RETURN
         END IF
      END IF
      GFLF => lst%sub(i)
      
      fName = GFLF%val
      fVal  = fileType(fName)

      io%d = TRIM(lst%ping(cmnd,GFLF))//
     2   " Opened file from path <"//TRIM(fName)//">"
      
      RETURN
      END FUNCTION GFLF

!####################################################################      
!     To test the implementations in this module and also to give you
!     a hint on how to use this class
      SUBROUTINE TEST_LSTMOD()

      LOGICAL l
      INTEGER i
      REAL(KIND=8) r, rV(3)
      CHARACTER(LEN=stdL) s
      TYPE(fileType) f
      TYPE(lstType) lst
      TYPE(lstType), POINTER :: lPtr

!     Create a dummy input file
      IF (cm%mas()) THEN
         i = 314
         OPEN (i,FILE='.dummy.mfs')
         WRITE(i,'(A)') "A logical input: F"
         WRITE(i,'(A)') "An integer input: 12"
         WRITE(i,'(A)') "A real input: 123D0"
         WRITE(i,'(A)') "A string input: Hello"
         WRITE(i,'(A)') "A vector input: 1D0 2D0 3D0"
         WRITE(i,'(A)') "A path to a file: ../../foo.dat"
         CLOSE(i)
      END IF

      lst = lstType(".dummy.mfs")
      io%o = "newLst: "//CLR("(PASSED)",3)
      lPtr => lst%get(l,"A logical input",1)
      IF (l) io%e = "Issue with logical read"
      io%o = "Logical read: "//CLR("(PASSED)",3)
      
      lPtr => lst%get(i,"An integer input",1)
      IF (i.NE.12) io%e = "Issue with integer read"
      io%o = "Integer read: "//CLR("(PASSED)",3)
      
      lPtr => lst%get(r,"A real input",1)
      IF (r.NE.123D0) io%e = "Issue with real read"
      io%o = "Real read: "//CLR("(PASSED)",3)
      
      lPtr => lst%get(s,"A string input",1)
      IF (s.NE."hello") io%e = "Issue with string read"
      io%o = "String read: "//CLR("(PASSED)",3)
      
      lPtr => lst%get(rV,"A vector input",1)
      IF (ANY(rV.NE.(/1D0,2D0,3D0/))) io%e = "Issue with vector read"
      io%o = "Vector read: "//CLR("(PASSED)",3)
      
      lPtr => lst%get(f,"A path to a file",1)
      IF (f%name().NE."../../foo.dat") io%e = "Issue with path read"
      io%o = "Path read: "//CLR("(PASSED)",3)

      IF (cm%mas()) THEN
         OPEN(1,FILE='.dummy.mfs')
         CLOSE(1,STATUS='DELETE')
      END IF

      RETURN
      END SUBROUTINE TEST_LSTMOD
      END MODULE LSTMOD
