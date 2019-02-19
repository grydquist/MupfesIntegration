!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Collection of all implemented equations. Also, implements routines
!     for reading and writing all variables. 
!      
!--------------------------------------------------------------------
      MODULE AEQMOD
      USE ADRMOD
      USE LEDMOD
      USE FSIMOD
      USE BBOMOD
      USE PRTMOD
      IMPLICIT NONE

!     Number of equations
      INTEGER nEq
!     All data related to equations are stored in this container
      TYPE(cEqType), POINTER :: eq(:)

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
      SUBROUTINE READ_AEQ(dmn,lst)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT), TARGET :: lst

      LOGICAL flag(dmn%nMsh)
      INTEGER iEq
      CHARACTER(LEN=stdL) ctmp
      TYPE(lstType), POINTER :: lPtr
      TYPE(insType), POINTER :: ns => NULL()
      TYPE(fsiType), POINTER :: fs => NULL()
      CLASS(eqType), POINTER :: eqPtr

!     Reading equations
      nEq = lst%srch("Add equation",1)
      io%o = " Number of equations: "//nEq
      ALLOCATE(eq(nEq))
      DO iEq=1, nEq
         lPtr => lst%get(ctmp,"Add equation",iEq)
         SELECT CASE (ctmp)
         CASE('adr')
            IF (ASSOCIATED(ns)) THEN
               ALLOCATE(eq(iEq)%s,SOURCE=adrType(dmn,lPtr,ns%U))
            ELSE
               ALLOCATE(eq(iEq)%s,SOURCE=adrType(dmn,lPtr))
            END IF
         CASE('led')
            ALLOCATE(eq(iEq)%s,SOURCE=ledType(dmn,lPtr))
         CASE('ins')
            ALLOCATE(eq(iEq)%s,SOURCE=insType(dmn,lPtr))
            eqPtr => eq(iEq)%s
            SELECT TYPE (eqPtr)
            TYPE IS (insType)
               ns => eqPtr
            END SELECT
         CASE('svk')
            ALLOCATE(eq(iEq)%s,SOURCE=svkType(dmn,lPtr))
         CASE('bbo')
            IF (.NOT.ASSOCIATED(ns)) io%e = "BBO must come after INS"
            ALLOCATE(eq(iEq)%s,SOURCE=bboType(ns, lPtr))
         CASE('fsi')
            ALLOCATE(eq(iEq)%s,SOURCE=fsiType(dmn, lPtr))
            eqPtr => eq(iEq)%s
            SELECT TYPE (eqPtr)
            TYPE IS (fsiType)
               fs => eqPtr
            END SELECT
         CASE('msh')
            IF (ASSOCIATED(fs)) THEN
               flag = fs%subPtr .EQ. 2
               ALLOCATE(eq(iEq)%s,SOURCE=ledType(dmn, lPtr,
     2            isMsh=.TRUE., g=fs%var(1), flag=flag))
            ELSE
               ALLOCATE(eq(iEq)%s,SOURCE=ledType(dmn, lPtr, 
     2            isMsh=.TRUE.))
            END IF
         CASE('prt')
            IF (.NOT.ASSOCIATED(ns)) io%e = "PRT must come after INS"
            ALLOCATE(eq(iEq)%s,SOURCE=prtType(dmn,lPtr,ns%U,ns%mat))
         CASE DEFAULT
            io%e = "Equation type "//TRIM(ctmp)//" is not defined"
         END SELECT
!     Initializing the remaining structures
         CALL eq(iEq)%s%ini()
      END DO
      IF (ASSOCIATED(fs) .AND. .NOT.ASSOCIATED(dmn%Um)) io%w =
     2   "FSI equation is solved without mesh motion (hope you know"//
     3   " what you're doing)!"

      RETURN
      END SUBROUTINE READ_AEQ
!---------------------------------------------------------------------
      SUBROUTINE RW_AEQ(dmn,f)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(fileType), INTENT(INOUT) :: f
      
      INTEGER iEq, i, pos 

      pos = 1
      IF (f%isRead()) THEN
         CALL cm%read(i,f%id(),pos)
         IF (i .NE. cm%np()) io%e = "Restart file: incompatible cm%np"
         CALL cm%read(i,f%id(),pos)
         IF (i .NE. nEq) io%e = "Restart file: incompatible nEq"
         CALL cm%read(i,f%id(),pos)
         IF (i .NE. dmn%dof) io%e = "Restart file: incompatible dmn%dof"
         CALL cm%read(cTS,f%id(),pos)
         CALL cm%read(time,f%id(),pos)
         DO iEq=1, nEq
            CALL cm%read(eq(iEq)%s%zNorm,f%id(),pos)
            IF (ALLOCATED(eq(iEq)%s%cbc%xn)) 
     2         CALL cm%read(eq(iEq)%s%cbc%xn,f%id(),pos)
            DO i=1, eq(iEq)%s%nVar
               CALL cm%read(eq(iEq)%s%var(i)%An%sP,f%id(),pos)
               CALL cm%read(eq(iEq)%s%var(i)% n%sP,f%id(),pos)
               CALL cm%read(eq(iEq)%s%var(i)%Dn%sP,f%id(),pos)
            END DO
         END DO
      ELSE
         CALL cm%write(cm%np(),f%id(),pos)
         CALL cm%write(nEq,f%id(),pos)
         CALL cm%write(dmn%dof,f%id(),pos)
         CALL cm%write(cTS,f%id(),pos)
         CALL cm%write(time,f%id(),pos)
         DO iEq=1, nEq
            CALL cm%write(eq(iEq)%s%zNorm,f%id(),pos)
            IF (ALLOCATED(eq(iEq)%s%cbc%xn)) 
     2         CALL cm%write(eq(iEq)%s%cbc%xn,f%id(),pos)
            DO i=1, eq(iEq)%s%nVar
               CALL cm%write(eq(iEq)%s%var(i)%An%sP,f%id(),pos)
               CALL cm%write(eq(iEq)%s%var(i)% n%sP,f%id(),pos)
               CALL cm%write(eq(iEq)%s%var(i)%Dn%sP,f%id(),pos)
            END DO
         END DO
      END IF
      CALL f%close()

      RETURN
      END SUBROUTINE RW_AEQ
      END MODULE AEQMOD
