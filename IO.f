!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This module is used for I/O. For usage see TEST_IOMOD at the end
!     of this file. 
!      
!---------------------------------------------------------------------

      MODULE IOMOD
      USE ISO_FORTRAN_ENV
      USE UTILMOD
      IMPLICIT NONE

!     Channel type, used in I/O
      TYPE, ABSTRACT :: chnlType
         PRIVATE
!        Whether it is open to the screen
         LOGICAL :: oTS = .FALSE.
!        Whether it is open to the file
         LOGICAL :: oTF = .FALSE.
!        File ID
         INTEGER fId
!        File address
         CHARACTER(LEN=stdL) :: fName = "histor"
!        Channel tag
         CHARACTER(LEN=stdL) :: tag = ""
      CONTAINS 
!        Closes the channel
         PROCEDURE :: close => closeChnl
!        To send a string to a channel
         PROCEDURE(printChnlIf), DEFERRED :: print
         GENERIC :: ASSIGNMENT(=) => print
!        To adjust oTF and oTS status
         PROCEDURE :: set => setChnl
!        TO be used by sub-classes. To create a new channel
         PROCEDURE, PRIVATE :: new => newChnl
!        To be used by sub-classes. Sends a message to channel.
         PROCEDURE, PRIVATE :: msg => msgChnl 
      END TYPE chnlType
 
      INTERFACE
         SUBROUTINE printChnlIf(chnl,msg)
            IMPORT :: chnlType
            CLASS(chnlType), INTENT(INOUT) :: chnl
            CHARACTER(LEN=*), INTENT(IN) :: msg
         END SUBROUTINE printChnlIf
      END INTERFACE
!---------------------------------------------------------------------
!     For standard input and output
      TYPE, EXTENDS(chnlType) :: stdType
      CONTAINS
         PROCEDURE :: print => printStd
         FINAL :: closeStd
      END TYPE stdType

      INTERFACE stdType
         PROCEDURE :: newStd
      END INTERFACE
!---------------------------------------------------------------------
!     For handling errors
      TYPE, EXTENDS(chnlType) :: errType
      CONTAINS
         PROCEDURE :: print => printErr
         FINAL :: closeErr
      END TYPE errType

      INTERFACE errType
         PROCEDURE :: newErr
      END INTERFACE
!---------------------------------------------------------------------
!     For handling debug messages
      TYPE, EXTENDS(chnlType) :: dbgType
      CONTAINS
         PROCEDURE :: print => printDbg
         FINAL :: closeDbg
      END TYPE dbgType

      INTERFACE dbgType
         PROCEDURE :: newDbg
      END INTERFACE
!---------------------------------------------------------------------
!     For handling warnings
      TYPE, EXTENDS(chnlType) :: wrnType
      CONTAINS
         PROCEDURE :: print => printWrn
!        Checking for exceptions
         PROCEDURE :: checkException => checkExceptionWrn
         FINAL :: closeWrn
      END TYPE wrnType

      INTERFACE wrnType
         PROCEDURE :: newWrn
      END INTERFACE
!---------------------------------------------------------------------
!     Grouping four channels into one variable
      TYPE ioType
!        Standard output
         TYPE(stdType) :: o
!        Error
         TYPE(errType) :: e
!        Debugging
         TYPE(dbgType) :: d
!        Warning
         TYPE(wrnType) :: w
!        Status file
         TYPE(stdType) :: s
      END TYPE ioType

      INTERFACE ioType
         PROCEDURE :: newIO
      END INTERFACE
!---------------------------------------------------------------------
!     Whether to use color in printing outputs
      LOGICAL :: pClr = .TRUE.
!     A general counter for file ID
      INTEGER, PRIVATE :: gFID = 314
!     Appended path to all files that are going to be saved
      CHARACTER(LEN=stdL) :: appPath = ""
!     A generic IO for my own use
      TYPE(ioType) :: io

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      
      SUBROUTINE newChnl(chnl)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(INOUT) :: chnl

      LOGICAL flag

      chnl%fName = TRIM(chnl%fName)//".dat"
      IF (chnl%tag .NE. "") chnl%fName = TRIM(chnl%fName)//"-"//chnl%tag

      gFID     = gFID + 1
      chnl%fId = gFID
      
      DO 
         INQUIRE(UNIT=chnl%fId, OPENED=flag)
         IF (.NOT.flag) EXIT
         chnl%fId = chnl%fId + 1
      END DO

      RETURN
      END SUBROUTINE newChnl
!---------------------------------------------------------------------
      SUBROUTINE setChnl(chnl,oTS,oTF)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(INOUT) :: chnl
      LOGICAL, INTENT(IN), OPTIONAL :: oTS, oTF

      IF (PRESENT(oTS)) chnl%oTS = oTS
      IF (PRESENT(oTF)) chnl%oTF = oTF

      RETURN
      END SUBROUTINE setChnl
!---------------------------------------------------------------------
      SUBROUTINE closeChnl(chnl)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(INOUT) :: chnl
      
      LOGICAL flag

      chnl%oTS   = .FALSE.
      chnl%oTF   = .FALSE.
      chnl%fName = "histor"
      chnl%tag   = ""
      INQUIRE(UNIT=chnl%fId, OPENED=flag)
      IF (flag) CLOSE(chnl%fId)

      RETURN
      END SUBROUTINE closeChnl
!---------------------------------------------------------------------
!     This routine is used for keeping track of what is printed on the
!     screen and history file
      SUBROUTINE msgChnl(chnl,msg)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: msg

      LOGICAL flag
      CHARACTER(LEN=stdL) sTmp, fName

      IF (.NOT.chnl%oTS .AND. .NOT.chnl%oTF) RETURN

!     Searching and removing all the color codes
      IF (chnl%oTF .OR. (chnl%oTS.AND.(.NOT.pClr))) sTmp = RMCLR(msg)

      IF (chnl%oTF) THEN
         INQUIRE(UNIT=chnl%fId, OPENED=flag)
         IF (.NOT.flag) THEN
            IF (appPath .NE. "") CALL SYSTEM("mkdir -p "//TRIM(appPath))
            fName = TRIM(appPath)//TRIM(chnl%fName)
            INQUIRE(FILE=fName, OPENED=flag)
            IF (.NOT.flag) THEN
               OPEN(chnl%fId, FILE=fName, POSITION="APPEND")
            ELSE
               INQUIRE(FILE=fName, NUMBER=chnl%fId)
            END IF
         END IF
         WRITE(chnl%fId,"(A)") TRIM(sTmp)
         CALL FLUSH(chnl%fId)
      END IF

      IF (chnl%oTS) THEN
         IF (chnl%tag .NE. "") THEN
            IF (pClr) THEN
               WRITE(*,"(A)") CLR(TRIM(chnl%tag)//" >>",4)//" "//
     2            TRIM(msg)
            ELSE
               WRITE(*,"(A)") TRIM(chnl%tag)//" >> "//TRIM(sTmp)
            END IF
         ELSE
            IF (pClr) THEN
               WRITE(*,"(A)") TRIM(msg)
            ELSE
               WRITE(*,"(A)") TRIM(sTmp)
            END IF
         END IF
         CALL FLUSH(OUTPUT_UNIT)
      END IF

      RETURN
      END SUBROUTINE msgChnl

!#####################################################################
!     Constructing a new std output
      FUNCTION newStd(fName,tag,oTS,oTF) RESULT(chnl)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fName, tag
      LOGICAL, INTENT(IN), OPTIONAL :: oTS, oTF
      TYPE(stdType) :: chnl

      IF (PRESENT(fName)) chnl%fName = fName
      IF (PRESENT(tag))   chnl%tag   = tag
      IF (PRESENT(oTS))   chnl%oTS   = oTS
      IF (PRESENT(oTF))   chnl%oTF   = oTF

      CALL chnl%new()

      RETURN
      END FUNCTION newStd
!---------------------------------------------------------------------
!     This routine is used for standard output handling
      SUBROUTINE printStd(chnl,msg)
      IMPLICIT NONE
      CLASS(stdType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: msg

      CALL chnl%msg(msg)

      END SUBROUTINE printStd
!---------------------------------------------------------------------
      SUBROUTINE closeStd(chnl)
      IMPLICIT NONE
      TYPE(stdType), INTENT(INOUT) :: chnl
      
      CALL chnl%close()

      RETURN
      END SUBROUTINE closeStd
!#####################################################################
!     Constructing a new error output
      FUNCTION newErr(fName,tag,oTS,oTF) RESULT(chnl)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fName, tag
      LOGICAL, INTENT(IN), OPTIONAL :: oTS, oTF
      TYPE(errType) :: chnl

      chnl%oTS = .TRUE.
      chnl%oTF = .TRUE.
      IF (PRESENT(fName)) chnl%fName = fName
      IF (PRESENT(tag))   chnl%tag   = tag
      IF (PRESENT(oTS))   chnl%oTS   = oTS
      IF (PRESENT(oTF))   chnl%oTF   = oTF

      CALL chnl%new()

      RETURN
      END FUNCTION newErr
!---------------------------------------------------------------------
!     This routine is used for error handling
      SUBROUTINE printErr(chnl,msg)
      IMPLICIT NONE
      CLASS(errType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: msg

      CHARACTER(LEN=LEN(msg)) sTmp

      IF (chnl%oTS .OR. chnl%oTF) THEN
         sTmp = RMCLR(msg)
         CALL chnl%msg("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//
     2      "!!!!!!!!")
         CALL chnl%msg("An ERROR occured. See below for more"//
     2      " detail")
         CALL chnl%msg(CLR("ERROR: "//TRIM(ADJUSTL(sTmp))))
         CALL chnl%msg("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//
     2      "!!!!!!!!")
      END IF

      ERROR STOP "All processors are forced to stop by a fatal error"
      END SUBROUTINE printErr
!---------------------------------------------------------------------
      SUBROUTINE closeErr(chnl)
      IMPLICIT NONE
      TYPE(errType), INTENT(INOUT) :: chnl
      
      CALL chnl%close()

      RETURN
      END SUBROUTINE closeErr
!#####################################################################
!     Constructing a new debug output
      FUNCTION newDbg(fName,tag,oTS,oTF) RESULT(chnl)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fName, tag
      LOGICAL, INTENT(IN), OPTIONAL :: oTS, oTF
      TYPE(dbgType) :: chnl

      IF (PRESENT(fName)) chnl%fName = fName
      IF (PRESENT(tag))   chnl%tag   = tag
      IF (PRESENT(oTS))   chnl%oTS   = oTS
      IF (PRESENT(oTF))   chnl%oTF   = oTF

      CALL chnl%new()

      RETURN
      END FUNCTION newDbg
!---------------------------------------------------------------------
!     This routine is used for handling debugging messages
      SUBROUTINE printDbg(chnl,msg)
      IMPLICIT NONE
      CLASS(dbgType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: msg

      IF (chnl%oTS .OR. chnl%oTF)
     2   CALL chnl%msg(CLR(" DEBUG:",3)//" "//ADJUSTL(msg))
      
      RETURN
      END SUBROUTINE printDbg
!---------------------------------------------------------------------
      SUBROUTINE closeDbg(chnl)
      IMPLICIT NONE
      TYPE(dbgType), INTENT(INOUT) :: chnl
      
      CALL chnl%close()

      RETURN
      END SUBROUTINE closeDbg
!#####################################################################
!     Constructing a new warning output
      FUNCTION newWrn(fName,tag,oTS,oTF) RESULT(chnl)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fName, tag
      LOGICAL, INTENT(IN), OPTIONAL :: oTS, oTF
      TYPE(wrnType) :: chnl

      chnl%oTS = .TRUE.
      chnl%oTF = .TRUE.
      IF (PRESENT(fName)) chnl%fName = fName
      IF (PRESENT(tag))   chnl%tag   = tag
      IF (PRESENT(oTS))   chnl%oTS   = oTS
      IF (PRESENT(oTF))   chnl%oTF   = oTF

      CALL chnl%new()

      RETURN
      END FUNCTION newWrn
!---------------------------------------------------------------------
!     This routine is used for warning messages
      SUBROUTINE printWrn(chnl,msg)
      IMPLICIT NONE
      CLASS(wrnType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: msg
      
      IF (chnl%oTS .OR. chnl%oTF) THEN
         IF (pClr) THEN
            CALL chnl%msg(CLR(" WARNING:",4)//" "//ADJUSTL(msg))
         ELSE
            CALL chnl%msg("!!! WARNING: "//ADJUSTL(msg))
         END IF
      END IF

      RETURN
      END SUBROUTINE printWrn
!---------------------------------------------------------------------
!     This is to find if any exception occures, commenting this out
!     untill fortran2003 is standard in all compilers
      SUBROUTINE checkExceptionWrn(wrn)
      USE IEEE_EXCEPTIONS 
      CLASS(wrnType), INTENT(INOUT) :: wrn

!     I am assuming nExCk is 5, if a compilation error occured, you need
!     to adjust the following arrays accordingely
      INTEGER, PARAMETER :: nExCk = SIZE(IEEE_ALL) 
      LOGICAL, PARAMETER :: check(nExCk) = 
     2   (/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE./)
      CHARACTER(LEN=stdL), PARAMETER :: ieWarn(nExCk) = 
     2   (/"Overflow", "Divide by zero", "Invalid arithmetic operation",
     3   "Underflow", "Inexact operation"/)
      LOGICAL, SAVE :: iniSet = .TRUE., sprtFlag(nExCk)
      LOGICAL fg
      INTEGER i
      REAL(KIND=8) r

      IF (iniSet) THEN
         iniSet = .FALSE.
         DO i=1, nExCk
            IF (.NOT.check(i)) CYCLE
            sprtFlag(i) = IEEE_SUPPORT_FLAG(IEEE_ALL(i), r)
         END DO
         
         fg = .FALSE.
         DO i=1, nExCk
            IF (.NOT.check(i)) CYCLE
            IF (sprtFlag(i)) CALL IEEE_SET_HALTING_MODE(IEEE_ALL(i), fg)
         END DO
      END IF

      DO i=1, nExCk
         IF (.NOT.check(i)) CYCLE
         IF (sprtFlag(i)) THEN
            CALL IEEE_GET_FLAG(IEEE_ALL(i), fg)
            IF (fg) THEN
               wrn = ieWarn(i)
               CALL IEEE_SET_FLAG(IEEE_ALL(i), .FALSE.)
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE checkExceptionWrn
!---------------------------------------------------------------------
      SUBROUTINE closeWrn(chnl)
      IMPLICIT NONE
      TYPE(wrnType), INTENT(INOUT) :: chnl
      
      CALL chnl%close()

      RETURN
      END SUBROUTINE closeWrn
!#####################################################################
      FUNCTION newIO(master, fName, tag) RESULT(io)
      IMPLICIT NONE
      LOGICAL, INTENT(IN), OPTIONAL :: master
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fName, tag
      TYPE(ioType) :: io
      
      LOGICAL flag

!     In case of multi processor, we only want the master to print stuff
!     into the screen/files
      flag = .TRUE.
      IF (PRESENT(master)) flag = master

      io%o = stdType(fName=fName,tag=tag,oTS=flag,oTF=flag)
      io%e = errType(fName=fName,tag=tag)
      io%w = wrnType(fName=fName,tag=tag)
      io%d = dbgType(fName=fName,tag=tag)
      io%s = stdType(fName=fName,tag=tag,oTF=flag)

      RETURN
      END FUNCTION newIO
!#####################################################################
!     To test the current implementation: Initializing an IO type, 
!     then doing things with it and closing it at the end
      SUBROUTINE TEST_IOMOD()
      
      io%o = "Standard channel:"
      io%o = "Hello world!"
      io%o = "Debugg channel:"
      io%d = "Don't expect to see this anywhere."
      io%o = "Warning channel:"
      io%w = "This is a warning message."
      io%o = "File output channel:"
      io%s = "This ought to show up only in the histor file."
      io%o = "Error channel:"
      io%o = "Huh! you don't want me to terminate the program"
      io%o = CLR("(PASSED)",3)

      RETURN
      END SUBROUTINE TEST_IOMOD

      END MODULE IOMOD

