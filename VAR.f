!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     Contains variables class that are constructed based on domain
!     class. Two types of variables are supported: scalar and vector.
!     For usage and testing this class, see TEST_VARMOD at the end of 
!     this file. 
!      
!---------------------------------------------------------------------
      MODULE VARMOD
      USE MSHMOD
      IMPLICIT NONE

      TYPE :: dmnType
!        Whether the mesh is initialized
         LOGICAL :: izd = .FALSE.
!        Number of meshes 
         INTEGER :: nMsh = 0
!        Number of faces couples projected
         INTEGER :: nPrj = 0
!        Global total number of nodes (nodes of all meshes)
         INTEGER :: gdof = 0
!        Total number of nodes (nodes of all meshes in this processor)
         INTEGER :: dof = 0
!        A varibale to keep track of projected faces. 
         CHARACTER(LEN=stdL), ALLOCATABLE :: prjFa(:,:)
!        Local to global dof  [dof --> gdof]
         INTEGER, ALLOCATABLE :: guP(:)
!        All the meshes are stored here
         TYPE(mshType), ALLOCATABLE :: msh(:)
!        To keep track of mesh motion
         TYPE(gVarType), POINTER, PUBLIC :: Um => NULL()
      CONTAINS
!        Finds the mesh and face index given the face name
         PROCEDURE :: find => findDmn
!        Initializes by partitioning all meshes belonging to this domain
         PROCEDURE :: ini => iniDmn
!        To deallocate an element
         PROCEDURE :: free => freeDmn
!        Does the initial setup ... 
         PROCEDURE, PRIVATE :: setup => setupDmn
!        Finds corresponding nodes to assign them to a single unknown
         PROCEDURE, PRIVATE :: setPrj => setPrjDmn 
      END TYPE dmnType

      INTERFACE dmnType
         PROCEDURE :: newDmn, newDmnFl
      END INTERFACE dmnType
!---------------------------------------------------------------------
      TYPE :: varType
!        The degrees of freedom per node
         INTEGER :: dof = 1
!        The name of the variable
         CHARACTER(LEN=stdL) :: name = DEFAULT_NAME
!        Pointer to the domain
         TYPE(dmnType), POINTER :: dmn => NULL()
!        Generic storage for variable
         REAL(KIND=8), POINTER :: sP(:) => NULL()
!        For scalar operations 
         REAL(KIND=8), POINTER :: s(:) => NULL()
!        For vector operations
         REAL(KIND=8), POINTER, CONTIGUOUS :: v(:,:) => NULL()
      CONTAINS
!        Deallocates variables
         PROCEDURE :: free => freeVar
      END TYPE varType

      INTERFACE varType
         PROCEDURE :: newVar
      END INTERFACE varType
!---------------------------------------------------------------------
!     For Generalized-alpha time integration and also tracking
!     integration and time derivative of a variable
!        The extended variable is the current version of the variable
      TYPE, EXTENDS(varType) :: gVarType
!        The old time derivative of the variable (e.g. acceleration)
         TYPE(varType) :: Ao
!        The cCurrent time derivative of the variable
         TYPE(varType) :: A
!        The new time derivative of the variable
         TYPE(varType) :: An
!        The older version of the variables (e.g. velocity)
         TYPE(varType) :: o
!        The new version of the variables
         TYPE(varType) :: n
!        The old version of the integrated variable (e.g. dissplacement)
         TYPE(varType) :: Do
!        The current version of the integrated variable
         TYPE(varType) :: D
!        The new version of the integrated variables
         TYPE(varType) :: Dn
!        Outside contributions to variable
         TYPE(varType) :: OC
      CONTAINS
!        Deallocating all variables
         PROCEDURE :: free => freeGVar
      END TYPE gVarType
      
      INTERFACE gVarType
         PROCEDURE :: newGVar
      END INTERFACE gVarType

      CONTAINS

!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newDmn(nMsh, nPrj) RESULT(dmn)
      INTEGER, INTENT(IN) :: nMsh, nPrj
      TYPE(dmnType) :: dmn

      ALLOCATE(dmn%msh(nMsh), dmn%prjFa(2,nPrj))
      dmn%nMsh  = nMsh
      dmn%nPrj  = nPrj
      dmn%prjFa = DEFAULT_NAME

      RETURN
      END FUNCTION newDmn
!---------------------------------------------------------------------
      FUNCTION newDmnFl(lst) RESULT(dmn)
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(dmnType) :: dmn

      INTEGER nMsh, nPrj, iM, iPrj
      CHARACTER(LEN=stdL) ctmp
      TYPE(lstType), POINTER :: lPtr, lPP

      nMsh = lst%srch("Add mesh",ll=1)
      io%o = " Number of meshes: "//nMsh
      nPrj = lst%srch("Add projection")
      IF (nPrj .NE. 0) io%o  = " Number of projection: "//nPrj
      dmn = newDmn(nMsh, nPrj)

      DO iM=1, nMsh
         lPtr => lst%get(ctmp,"Add mesh",iM)
         io%o = " Reading mesh <"//CLR(TRIM(ctmp),3)//">"
         dmn%msh(iM) = mshType(lPtr, ctmp)
      END DO

      DO iPrj=1, nPrj
         lPP  => lst%get(dmn%prjFa(1,iPrj),"Add projection",iPrj)
         lPtr => lPP%get(dmn%prjFa(2,iPrj),"Project from face",1)
      END DO

      RETURN
      END FUNCTION newDmnFl
!---------------------------------------------------------------------
      SUBROUTINE setupDmn(dmn)
      CLASS(dmnType), INTENT(INOUT), TARGET :: dmn

      INTEGER iM, iFa, jM, jFa, i, a, Ac
      CHARACTER(LEN=stdL) ctmp
      TYPE(stackType) avNds
      
!     Making sure all meshes are assigned and face names are unique. 
!     To be call sequentially
      DO iM=1, dmn%nMsh
         IF (dmn%msh(iM)%nNo .EQ. 0) io%e = "All meshes must be"
     2      //"assigned to the domain before a call to setup"
         ALLOCATE(dmn%msh(iM)%uP(dmn%msh(iM)%nNo))
         dmn%msh(iM)%uP = 0
         DO iFa=1, dmn%msh(iM)%nFa
            ctmp = dmn%msh(iM)%fa(iFa)%name
            DO jM=1, iM
               DO jFa=1, dmn%msh(jM)%nFa
                  IF (ctmp .EQ. dmn%msh(jM)%fa(jFa)%name .AND. 
     2               (jM.NE.iM .OR. jFa.NE.iFa)) THEN
                     io%e = "Face names should be unique"
                  END IF
               END DO
            END DO
         END DO
      END DO

!     Examining the existance of projection faces and setting %uP.
!     Seting gdof and recounting nodes that are not duplicated
      dmn%gdof = 0
      CALL dmn%setPrj(avNds)
      DO iM=1, dmn%nMsh
         DO a=1, dmn%msh(iM)%nNo
            IF (dmn%msh(iM)%uP(a) .EQ. 0) THEN
               IF (avNds%pull(i)) THEN
                  dmn%msh(iM)%uP(a) = i
               ELSE
                  dmn%gdof          = dmn%gdof + 1
                  dmn%msh(iM)%uP(a) = dmn%gdof
               END IF
            END IF
         END DO
      END DO
      IF (avNds%pull(i)) io%e = "There are still nodes in the stack"
      CALL avNds%free()

!     Constructing face uP
      DO iM=1, dmn%nMsh
         DO iFa=1, dmn%msh(iM)%nFa
            ALLOCATE(dmn%msh(iM)%fa(iFa)%uP(dmn%msh(iM)%fa(iFa)%nNo))
            dmn%msh(iM)%fa(iFa)%uP = 0
            DO a=1, dmn%msh(iM)%fa(iFa)%nNo
               Ac = dmn%msh(iM)%fa(iFa)%gN(a)
               dmn%msh(iM)%fa(iFa)%uP(a) = dmn%msh(iM)%uP(Ac)
            END DO
         END DO
      END DO

      io%o = CLR(" Mesh data imported successfully",3)

      RETURN
      END SUBROUTINE setupDmn
!---------------------------------------------------------------------
!     This routines associates two faces with each other and sets gN
      SUBROUTINE setPrjDmn(dmn,avNds)
      CLASS(dmnType), INTENT(INOUT) :: dmn
      TYPE(stackType), INTENT(OUT) :: avNds
      
      INTEGER i, j, k, iM, jM, kM, iFa, jFa, ia, ja, nPrj, iPrj, nStk
      TYPE(stackType) lPrj
      TYPE(stackType), ALLOCATABLE :: stk(:)

!     Here calculating an upper limit for the number required stacks
      nStk = 0
      nPrj = dmn%nPrj
      DO iPrj=1, nPrj
         CALL dmn%find(dmn%prjFa(1,iPrj), iM, iFa)
         nStk = nStk + dmn%msh(iM)%fa(iFa)%nNo
      END DO
      ALLOCATE(stk(nStk))

      DO iPrj=1, nPrj
         CALL dmn%find(dmn%prjFa(1,iPrj), iM, iFa)
         CALL dmn%find(dmn%prjFa(2,iPrj), jM, jFa)
         lPrj = dmn%msh(iM)%matchFace(dmn%msh(jM), iFa, jFa)
         DO
!     First in, last out
            IF (.NOT.lPrj%pull(ja)) EXIT
            IF (.NOT.lPrj%pull(ia)) EXIT
            i = dmn%msh(iM)%uP(ia)
            j = dmn%msh(jM)%uP(ja)
            IF (i .EQ. 0) THEN
               IF (j .EQ. 0) THEN
!     Since neither of them have value, I will add a new node and both
!     of them to the stack
                  IF (.NOT.avNds%pull(k)) THEN
                     dmn%gdof = dmn%gdof + 1
                     k        = dmn%gdof
                  END IF
                  dmn%msh(iM)%uP(ia) = k
                  dmn%msh(jM)%uP(ja) = k
                  CALL stk(k)%push((/iM,ia,jM,ja/))
               ELSE
!     This is the case one of them has already been assigned. So just
!     using that value for the other one
                  dmn%msh(iM)%uP(ia) = j
                  CALL stk(j)%push((/iM,ia/))
               END IF
            ELSE
               IF (j .EQ. 0) THEN
                  dmn%msh(jM)%uP(ja) = i
                  CALL stk(i)%push((/jM,ja/))
               ELSE
!     Since they are both already have assigned values, I will move the
!     nodes from stack with bigger ID, j, to the other stack, i.
                  IF (i .EQ. j) CYCLE
                  IF (i .GT. j) THEN
                     k = i
                     i = j
                     j = k
                  END IF
                  DO
                     IF (.NOT.stk(j)%pull(ja)) EXIT
                     IF (.NOT.stk(j)%pull(kM)) EXIT
                     dmn%msh(kM)%uP(ja) = i
                     CALL stk(i)%push((/kM,ja/))
                  END DO
                  CALL avNds%push(j)
               END IF
            END IF
         END DO
!     Since all nodes are added, I remove all the members
         CALL lPrj%free()
      END DO

      RETURN
      END SUBROUTINE setPrjDmn
!---------------------------------------------------------------------
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE findDmn(dmn, name, iM, iFa)
      CLASS(dmnType), INTENT(IN) :: dmn
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(OUT) :: iM, iFa

      CHARACTER(LEN=LEN(TRIM(ADJUSTL(name)))) fName
      
      fName = LOWER(TRIM(ADJUSTL(name)))
      iFa   = 0
      MY_LOOP : DO iM=1, dmn%nMsh
         DO iFa=1, dmn%msh(iM)%nFa
            IF (LOWER(TRIM(ADJUSTL(dmn%msh(iM)%fa(iFa)%name))) .EQ. 
     2         fName) EXIT MY_LOOP
         END DO
      END DO MY_LOOP

      IF (iM .GT. dmn%nMsh) io%e = "Unable to find face <"//fName//">"

      RETURN
      END SUBROUTINE findDmn
!---------------------------------------------------------------------
      SUBROUTINE iniDmn(dmn)
      CLASS(dmnType), INTENT(INOUT) :: dmn

      INTEGER nMsh, iM, i, iFa
      INTEGER, ALLOCATABLE :: gmtl(:)
      REAL(KIND=8), ALLOCATABLE :: iWgt(:), wgt(:,:), wrk(:)
 
      IF (dmn%izd) THEN
         io%w = "Domain already initialzied"
         RETURN
      END IF
      dmn%izd = .TRUE.

      IF (cm%mas()) CALL dmn%setup()

!     Constructing data structures one by one
      CALL cm%bcast(dmn%gdof)

!     wgt and wrk are the assigned portion of each mesh to the each 
!     processor.
      nMsh = dmn%nMsh
      ALLOCATE(wgt(nMsh,cm%np()), wrk(nMsh), iWgt(cm%np()))

!     Here is rough estimation of how each mesh should be splited
!     between processors
      wrk = REAL(dmn%msh%nNo,8)/REAL(dmn%gdof,8)
      CALL cm%bcast(wrk)
      wgt = SPLIT_JOBS(nMsh, cm%np(), wrk)

!     First partitioning the meshes 
      IF (cm%seq()) THEN
!     Assigning all DOF fields
         dmn%dof = dmn%gdof
         ALLOCATE(dmn%guP(dmn%dof), gmtl(dmn%gdof))
         dmn%guP = (/(i, i=1, dmn%dof)/)
         gmtl    = (/(i, i=1, dmn%gdof)/)
      ELSE
         ALLOCATE(gmtl(dmn%gdof))
         gmtl = 0
      END IF

      DO iM=1, nMsh
         io%d = "Partitioning mesh "//iM
         iWgt = wgt(iM,:)/SUM(wgt(iM,:))
         CALL dmn%msh(iM)%ini(iWgt, dmn%dof, dmn%guP, gmtl)
      END DO

!     Partitioning all faces
      DO iM=1, nMsh
         DO iFa=1, dmn%msh(iM)%nFa
            CALL dmn%msh(iM)%partFa(iFa, gmtl)
         END DO
         DEALLOCATE(dmn%msh(iM)%otn, dmn%msh(iM)%gtl)
      END DO
      io%o = " Domain is partitioned successfully."

      RETURN
      END SUBROUTINE iniDmn
!---------------------------------------------------------------------
      SUBROUTINE freeDmn(this)
      CLASS(dmnType) :: this

      INTEGER iM

      IF (ALLOCATED(this%prjFa)) DEALLOCATE(this%prjFa)
      IF (ALLOCATED(this%guP))   DEALLOCATE(this%guP)
      IF (ALLOCATED(this%msh)) THEN
         DO iM=1, this%nMsh
            CALL this%msh(iM)%free()
         END DO
         DEALLOCATE(this%msh)
      END IF

      this%izd  = .FALSE.
      this%nMsh = 0
      this%nPrj = 0
      this%gdof = 0
      this%dof  = 0
     
      RETURN 
      END SUBROUTINE freeDmn

!#####################################################################
      
      FUNCTION newVar(dof, name, dmn) RESULT(var)
      INTEGER, INTENT(IN) :: dof
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      TYPE(varType) :: var

      IF (.NOT.dmn%izd) io%e = "newVar: Un-initialized domain"

      var%dof  = dof
      var%dmn => dmn
      var%name = ADJUSTL(name)
      
      ALLOCATE(var%sP(dof*dmn%dof))
      var%sP = 0D0

      IF (dof .EQ. 1) THEN
         var%s => var%sP
      ELSE IF (dof .EQ. nsd) THEN
         var%v(1:nsd,1:dmn%dof) => var%sP
      ELSE
         io%e = "newVar: only defined for dof=1 and dof=nsd"
      END IF

      RETURN
      END FUNCTION newVar
!---------------------------------------------------------------------
      SUBROUTINE freeVar(var)
      CLASS(varType), INTENT(INOUT) :: var

      var%dof  = 1
      var%name = DEFAULT_NAME
      var%dmn => NULL()
      var%s   => NULL()
      var%v   => NULL()

      IF (ASSOCIATED(var%sP)) DEALLOCATE(var%sP)
      
      RETURN
      END SUBROUTINE freeVar
!#####################################################################
      FUNCTION newGVar(dof, name, dmn) RESULT(var)
      INTEGER, INTENT(IN) :: dof
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      TYPE(gVarType) :: var

      var%Ao = varType(dof, name, dmn)
      var%A  = varType(dof, name, dmn)
      var%An = varType(dof, name, dmn)
      var%o = varType(dof, name, dmn)
      var%varType = varType(dof, name, dmn)
      var%n = varType(dof, name, dmn)
      var%Do = varType(dof, name, dmn)
      var%D  = varType(dof, name, dmn)
      var%Dn = varType(dof, name, dmn)
      var%OC = varType(dof, name, dmn)

      RETURN
      END FUNCTION newGVar
!---------------------------------------------------------------------
      SUBROUTINE freeGVar(var)
      CLASS(gVarType), INTENT(INOUT) :: var

      CALL var%varType%free()

      CALL var%Ao%free()
      CALL var%A %free()
      CALL var%An%free()
      CALL var%o%free()
      CALL var%varType%free()
      CALL var%n%free()
      CALL var%Do%free()
      CALL var%D %free()
      CALL var%Dn%free()
      CALL var%OC%free()

      RETURN
      END SUBROUTINE freeGVar

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_VARMOD(nMsh, dmn, prs, vel, gPrs)
      INTEGER, INTENT(IN) :: nMsh
      TYPE(dmnType), INTENT(OUT), TARGET :: dmn
      TYPE(varType), INTENT(OUT) :: prs, vel
      TYPE(gVarType), INTENT(OUT) :: gPrs

!     creating a dummy mesh
      IF (cm%mas()) CALL CREATE_DUMMY_MESH('.dummy') 
!     num_of_msh, num_of_proj
      dmn = dmnType(nMsh,nMsh-1) 
      io%o = "newDmn: "//CLR("(PASSED)",3)
!     mesh_name, num_of_face
      dmn%msh(1) = mshType("Mesh_1",2) 
!     Reading mesh from vtk file
      CALL dmn%msh(1)%read(".dummy.vtk") 
!     Reading face from vtk
      CALL dmn%msh(1)%readFace(1,"Face_p1",".dummypFace.vtk") 
      CALL dmn%msh(1)%readFace(2,"Face_1",".dummyFace.vtk") 
      IF (nMsh .EQ. 2) THEN
         dmn%msh(2) = mshType("Mesh_2",2)
         CALL dmn%msh(2)%read(".dummy.vtk")
         CALL dmn%msh(2)%readFace(1,"Face_p2",".dummypFace.vtk") 
         CALL dmn%msh(2)%readFace(2,"Face_2",".dummynFace.vtk") 
         dmn%prjFa(:,1) = (/"Face_p1","Face_p2"/)
      END IF
!     Since I am done reading meshes ...
      CALL dmn%ini() 
      IF (nMsh.EQ.2 .AND. dmn%gdof.NE.12) io%e = "Issue with ini"
      IF (nMsh.EQ.1 .AND. dmn%gdof.NE.8) io%e = "Issue with ini"
      io%o = "dmn%ini: "//CLR("(PASSED)",3)
      IF (cm%slv().AND.ANY(dmn%guP.NE.(/2,3,10,4,12,9,1,11/)) .AND.
     1   nMsh.EQ.2) io%e = "Issue with part"
      io%o = "dmn%part: "//CLR("(PASSED)",3)
!     Testing varTypes
      prs = varType(1,'pressure',dmn)
      vel = varType(nsd,'velocity',dmn)
      IF (.NOT.ASSOCIATED(prs%s) .OR. ASSOCIATED(prs%v) .OR. 
     1    .NOT.ASSOCIATED(vel%v) .OR. ASSOCIATED(vel%s))
     2   io%e = "Issue with newVar"
      io%o = "newVar: "//CLR("(PASSED)",3)
      gPrs = gVarType(1,'pressure',dmn)
      IF (.NOT.ASSOCIATED(gPrs%s) .OR. ASSOCIATED(gPrs%v)) 
     1   io%e = "Issue with newGVar"
      io%o = "newGVar: "//CLR("(PASSED)",3)
      
      RETURN
      END SUBROUTINE TEST_VARMOD
      END MODULE VARMOD
