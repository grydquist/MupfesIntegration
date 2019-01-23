!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Creating the data structure and assembling LHS sparse matrix. Also
!     handles the calls to memLS and stores related variables. 
!     For usage and testing this class, see TEST_LSMOD at the end of
!     this file. 
!      
!--------------------------------------------------------------------
      MODULE LSMOD
      USE VARMOD
      IMPLICIT NONE
      INCLUDE "memLS.h"

      TYPE, EXTENDS(memLS_lsType) :: lsType
c         PRIVATE
!        Whether it is initialized
         LOGICAL :: izd = .FALSE.
!        Degrees of freedom
         INTEGER :: dof = 0
!        Number of unknowns per dof
         INTEGER, PUBLIC :: nU = 0
!        Number of nonzeros in this proc.
         INTEGER :: nnz = 0
!        Number of faces to be setup with LS
         INTEGER, PUBLIC :: nFa = 0
!        Column pointer (for sparse LHS matrix structure)
         INTEGER, ALLOCATABLE :: colPtr(:)
!        Row pointer (for sparse LHS matrix structure)
         INTEGER, ALLOCATABLE :: rowPtr(:)
!        RHS vector
         REAL(KIND=8), ALLOCATABLE :: R(:,:)
!        LHS matrix 
         REAL(KIND=8), ALLOCATABLE :: Val(:,:)
!        memLS data structure to produce LHS sparse matrix
         TYPE(memLS_lhsType) lhs
      CONTAINS
!        Adds a BC
         PROCEDURE :: addBc => addBcLs
!        Assembling a local LHS/RHS into the system
         PROCEDURE :: assem => assemLs
!        To initialize
         PROCEDURE :: ini => iniLs
!        Solve the linear system
         PROCEDURE :: solve => solveLs
!        Rsets R and Val to zero
         PROCEDURE :: reset => resetLs
!        Returns the solution
         PROCEDURE :: soln => solnLs
!        Returns izd status
         PROCEDURE :: isIzd => isIzdLs
!        Deallocates all variables
         PROCEDURE :: free => freeLs
!        Construct LHS sparse structure
         PROCEDURE, PRIVATE :: cLhs => cLhsLs
      END TYPE lsType

      INTERFACE lsType
         PROCEDURE :: newLs, newLsFl
      END INTERFACE lsType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
!     Creates a new ls instance provided linear solver algorithm and
!     number of unknowns per degrees of freedom
      FUNCTION newLs(lsAl, nFa, dmn) RESULT(ls)
      TYPE(lsType) :: ls
      INTEGER, INTENT(IN) :: lsAl, nFa
      TYPE(dmnType), INTENT(IN) :: dmn

      INTEGER gnnz
      TYPE(memLS_commuType) communicator
      
      IF (.NOT.dmn%izd) io%e = "newLs: Un-initialized domain"

      ls%nFa = nFa
      CALL memLS_LS_CREATE(ls, lsAl)
      
      io%o = " Constructing stiffness matrix sparse structure"
      CALL ls%cLhs(dmn)
       
      gnnz = cm%reduce(ls%nnz)
      io%o = " Total number of non-zeros in the LHS matrix: "//gnnz

      io%d = "Calling memLS_COMMU_CREATE"
      CALL memLS_COMMU_CREATE(communicator, cm%com())
      
      io%d = "Calling memLS_LHS_CREATE"
      CALL memLS_LHS_CREATE(ls%lhs, communicator, dmn%gdof, ls%dof, 
     2   ls%nnz, dmn%guP, ls%rowPtr, ls%colPtr, ls%nFa)

      RETURN
      END FUNCTION newLs
!---------------------------------------------------------------------
!     Creates a new ls instance from lst lst, using lsAl as default
!     linear solver algorithm.
      FUNCTION newLsFl(lst, lsAl, nFa, dmn) RESULT(ls)
      TYPE(lstType), INTENT(INOUT) :: lst
      INTEGER, INTENT(IN) :: lsAl, nFa
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lsType) :: ls

      INTEGER flsAl
      CHARACTER(LEN=stdL) ctmp
      TYPE(lstType), POINTER :: lPtr, lPL

      flsAl = lsAl
      lPL => lst%get(ctmp,"LS type")
      IF (ASSOCIATED(lPL)) THEN
         SELECT CASE(TRIM(ctmp)) 
         CASE('NS')
            flsAl = LS_TYPE_NS
         CASE('GMRES')
            flsAl = LS_TYPE_GMRES
         CASE('CG')
            flsAl = LS_TYPE_CG
         CASE DEFAULT
            io%e = TRIM(lst%ping("LS type",lPL))//" Undefined type"
         END SELECT
      END IF
      ls = newLs(flsAl, nFa, dmn)

      IF (ASSOCIATED(lPL)) THEN
         lPtr => lPL%get(ls%RI%mItr,"Max iterations",ll=1)
         lPtr => lPL%get(ls%GM%mItr,"NS-GM max iterations",ll=1)
         lPtr => lPL%get(ls%CG%mItr,"NS-CG max iterations",ll=1)
         
         lPtr => lPL%get(ls%RI%relTol,"Tolerance",lb=0D0,ub=1D0)
         lPtr => lPL%get(ls%GM%relTol,"NS-GM tolerance",lb=0D0,ub=1D0)
         lPtr => lPL%get(ls%CG%relTol,"NS-CG tolerance",lb=0D0,ub=1D0)
         
         lPtr =>lPL%get(ls%RI%absTol,"Absolute tolerance",lb=0D0,ub=1D0)
         ls%GM%absTol = ls%RI%absTol
         ls%CG%absTol = ls%RI%absTol

         lPtr => lPL%get(ls%RI%sD,"Krylov space dimension",ll=1)
         IF (ASSOCIATED(lPtr)) ls%GM%sD = ls%RI%sD
      END IF

      RETURN
      END FUNCTION newLsFl
!---------------------------------------------------------------------
!     To construct LHS matrix sparse structure
      SUBROUTINE cLhsLs(ls, dmn)
      CLASS(lsType), INTENT(INOUT) :: ls
      TYPE(dmnType), INTENT(IN) :: dmn

      LOGICAL flag
      INTEGER i, j, a, b, e, mnnzeic, rowN, colN, iM
      INTEGER, ALLOCATABLE :: uInd(:,:)

      flag = .FALSE.
      mnnzeic = 10*MAXVAL(dmn%msh%eNoN)
 003  mnnzeic = mnnzeic + MAX(5,mnnzeic/5)

      IF (ALLOCATED(uInd)) DEALLOCATE(uInd)
      ALLOCATE (uInd(mnnzeic,dmn%dof))
      uInd = 0
      DO iM=1, dmn%nMsh
         DO e=1, dmn%msh(iM)%nEl
            DO a=1, dmn%msh(iM)%eNoN
               rowN = dmn%msh(iM)%IEN(a,e)
               rowN = dmn%msh(iM)%uP(rowN)
               DO b=1, dmn%msh(iM)%eNoN
                  colN = dmn%msh(iM)%IEN(b,e)
                  colN = dmn%msh(iM)%uP(colN)
                  DO i=1, mnnzeic
!     If current entry is zero, then  fill it up
                     IF (uInd(i,rowN) .EQ. 0) THEN
                        uInd(i,rowN) = colN
                        EXIT
                     END IF
!     If current entry is still smaller, keep going
                     IF (colN .GT. uInd(i,rowN)) CYCLE
!     If there is already a similar entry, no point to add it again
                     IF (colN .EQ. uInd(i,rowN)) EXIT
!     If we are this point, then the current entry is bigger. Hence
!     we need to shift all the entries from here to the end of the list.
!     If list is full, we request a larger list, otherwise we
!     shift and add the item at the current entry position.
                     IF (uInd(mnnzeic,rowN) .NE. 0) GOTO 003
                     DO j=mnnzeic, i+1, -1
                        uInd(j,rowN) = uInd(j-1,rowN)
                     END DO
                     uInd(i,rowN) = colN
                     EXIT
                  END DO
!     In the following case the list is too short.
                  IF (i .GT. mnnzeic) THEN
                     flag = .TRUE.
                     GOTO 003
                  END IF
               END DO
            END DO
         END DO
      END DO
      IF (flag) io%o = "max(nnz) is increased to: "//mnnzeic

!     Finding number of non-zeros in colPtr vector
      ls%nnz = 0
      DO rowN=1, dmn%dof
         IF (uInd(1,rowN) .EQ. 0) THEN
            io%e = "Node "//rowN//" is isolated"
         END IF
         DO i = 1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               ls%nnz = ls%nnz + 1
            END IF
         END DO  
      END DO

      ls%dof = dmn%dof
!     Now constructing compact form of rowPtr and colPtr
      ALLOCATE(ls%colPtr(ls%nnz), ls%rowPtr(ls%dof+1))

      j  = 1
      ls%rowPtr(1) = 1
      DO rowN=1, ls%dof
         DO i=1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               ls%colPtr(j) = uInd(i,rowN)
               j = j + 1
            END IF
         END DO
         ls%rowPtr(rowN+1) = j
      END DO
      DEALLOCATE (uInd)

      RETURN
      END SUBROUTINE cLhsLs
!---------------------------------------------------------------------
!     Adds a BC to LS structure
      SUBROUTINE addBcLs(ls, lsPtr, nNo, bcGrp, gNodes, sVl)
      CLASS(lsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: lsPtr, nNo, bcGrp
      INTEGER, INTENT(IN) :: gNodes(nNo)
      REAL(KIND=8), INTENT(IN) :: sVl(nsd,nNo)
 
      CALL memLS_BC_CREATE(ls%lhs, lsPtr, nNo, nsd, bcGrp, gNodes, sVl)

      RETURN
      END SUBROUTINE addBcLs
!---------------------------------------------------------------------
!     Adds a BC to LS structure
      SUBROUTINE iniLs(ls, nU)
      CLASS(lsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: nU
      
      ls%izd = .TRUE.
      ls%nU  = nU
      ALLOCATE(ls%R(ls%nU,ls%dof), ls%Val(ls%nU*ls%nU,ls%nnz))
      
      RETURN
      END SUBROUTINE iniLs
!---------------------------------------------------------------------
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      PURE SUBROUTINE assemLs(ls, d, eqN, lK, lR)
      CLASS(lsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: d, eqN(d)
      REAL(KIND=8), INTENT(IN) :: lK(:,:,:), lR(:,:)

      INTEGER a, b, ptr, rowN, colN, left, right, n, i, j

      n = SIZE(lR,1)
      IF (n .EQ. ls%nU) THEN
         DO a=1, d
            rowN = eqN(a)
            ls%R(:,rowN) = ls%R(:,rowN) + lR(:,a)
            DO b=1, d
               colN = eqN(b)
               left  = ls%rowPtr(rowN)
               right = ls%rowPtr(rowN+1)
               ptr   = (right + left)/2
               DO WHILE (colN .NE. ls%colPtr(ptr))
                  IF (colN .GT. ls%colPtr(ptr)) THEN
                     left  = ptr
                  ELSE
                     right = ptr
                  END IF
                  ptr = (right + left)/2
               END DO
               ls%Val(:,ptr) = ls%Val(:,ptr) + lK(:,a,b)
            END DO
         END DO
      ELSE
         DO a=1, d
            rowN = eqN(a)
            ls%R(1:n,rowN) = ls%R(1:n,rowN) + lR(:,a)
            DO b=1, d
               colN = eqN(b)
               left  = ls%rowPtr(rowN)
               right = ls%rowPtr(rowN+1)
               ptr   = (right + left)/2
               DO WHILE (colN .NE. ls%colPtr(ptr))
                  IF (colN .GT. ls%colPtr(ptr)) THEN
                     left  = ptr
                  ELSE
                     right = ptr
                  END IF
                  ptr = (right + left)/2
               END DO
               DO i=1, n
                  j = (i-1)*ls%nU
                  ls%Val(j+1:j+n,ptr) = ls%Val(j+1:j+n,ptr) 
     2                                + lK((i-1)*n+1:i*n,a,b)
               END DO
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE assemLs
!---------------------------------------------------------------------
!     Solving the linear system of equations
      SUBROUTINE solveLs(ls, incL, res, isS)
      CLASS(lsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN), OPTIONAL :: incL(ls%nFa)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: res(ls%nFa)
      LOGICAL, INTENT(IN), OPTIONAL :: isS(ls%dof)
 
      IF (.NOT.ls%izd) io%e = "solveLs: Not initialized"
      io%d = "Solving linear system"
      io%d = "Get rid of isS, incL, and res"
      CALL memLS_SOLVE(ls%lhs, ls, ls%nU, ls%R, ls%Val, isS, incL, 
     2   res)

      RETURN
      END SUBROUTINE solveLs
!---------------------------------------------------------------------
      SUBROUTINE resetLs(ls)
      CLASS(lsType), INTENT(INOUT) :: ls

      IF (.NOT.ls%izd) io%e = "resetLs: Not initialized"
      ls%R   = 0D0
      ls%Val = 0D0

      RETURN
      END SUBROUTINE resetLs
!---------------------------------------------------------------------
      FUNCTION solnLs(ls) RESULT(R)
      CLASS(lsType), INTENT(IN) :: ls
      REAL(KIND=8), POINTER :: R(:)

      CALL PNT_TO(R, ls%R)

      RETURN
      END FUNCTION solnLs
!---------------------------------------------------------------------
      FUNCTION isIzdLs(ls) RESULT(izd)
      CLASS(lsType), INTENT(IN) :: ls
      LOGICAL :: izd

      izd = ls%izd

      RETURN
      END FUNCTION isIzdLs
!---------------------------------------------------------------------
!     Deallocates all structures
      SUBROUTINE freeLs(ls)
      CLASS(lsType), INTENT(INOUT) :: ls

!     Deallocating sparse matrix structures
      IF(ls%lhs%foc) CALL memLS_LHS_FREE(ls%lhs)
      
      IF (ALLOCATED(ls%colPtr)) DEALLOCATE(ls%colPtr)
      IF (ALLOCATED(ls%rowPtr)) DEALLOCATE(ls%rowPtr)
      IF (ALLOCATED(ls%R))      DEALLOCATE(ls%R)
      IF (ALLOCATED(ls%Val))    DEALLOCATE(ls%Val)

      ls%izd = .FALSE.
      ls%dof = 0
      ls%nU  = 0
      ls%nnz = 0
      ls%nFa = 0
      
      RETURN
      END SUBROUTINE freeLs
!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_LSMOD(dmn)
      TYPE(dmnType), INTENT(IN) :: dmn

      INTEGER, PARAMETER :: d = 3
      INTEGER eqN(d)
      REAL(KIND=8) lR(1,d), lK(1,d,d)
      TYPE(lsType) ls

      ls   = lsType(LS_TYPE_CG,1,dmn) ! .,nFa,.
      io%o = "newLs: "//CLR("(PASSED)",3)
      CALL ls%ini(1) ! number of unknowns
      io%o = "ls%ini: "//CLR("(PASSED)",3)
      IF (.NOT.ls%isIzd()) io%e = "Issue with ls%izd"
      io%o = "ls%isIzd: "//CLR("(PASSED)",3)
      IF (dmn%msh(1)%nNo .GE. d) THEN
         eqN = dmn%msh(1)%uP(1:d)
      ELSE IF (dmn%nMsh .GT. 1) THEN
         IF (dmn%msh(2)%nNo .GE. d) eqN = dmn%msh(2)%uP(1:d)
      ELSE
         io%e = "Unable to find a mesh with more "//d//" nodes"
      END IF
      lR  = 1D0
      lK  = 1D0
      CALL ls%assem(d,eqN,lK,lR)
      io%o = "ls%assem: "//CLR("(PASSED)",3)

!     Manually dealloctaing R and Val so that ls can be used later
      DEALLOCATE(ls%R, ls%val)
      ls%izd = .FALSE.

      RETURN
      END SUBROUTINE TEST_LSMOD
      END MODULE LSMOD
