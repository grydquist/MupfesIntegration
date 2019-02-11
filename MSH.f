!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Contains msh, i.e. mesh, class. To store, read, and partition
!     the mesh. You may also use the class procedures to map
!     local-to-processor to global representation of variables. For
!     usage and testing this class, see TEST_MSHMOD at the end of this
!     file. 
!      
!--------------------------------------------------------------------
      MODULE MSHMOD
      USE PMSHMOD
      IMPLICIT NONE

!     This is the container for a mesh
      TYPE, EXTENDS(pMshType) :: mshType
c         PRIVATE
!        Whether the mesh is initialized
         LOGICAL :: izd = .FALSE.
!        Global number of elements
         INTEGER :: gnEl
!        Global number of nodes
         INTEGER :: gnNo
!        Global node number [nNo --> gnNo]
         INTEGER, ALLOCATABLE :: gN(:)
!        Global to local node number (temporary) [gnNo --> nNo]
         INTEGER, ALLOCATABLE :: gtl(:)
!        Old to new element number (temporary) [gnEl --> gnEl]
         INTEGER, ALLOCATABLE :: otn(:)
!        Global connectivity array [eNoN x gnEl --> gnNo]
         INTEGER, ALLOCATABLE :: gIEN(:,:)
!        Distributuon elements between procs
         INTEGER, ALLOCATABLE :: eD(:)
!        Node to unknown pointer [nNo --> dof]
         INTEGER, ALLOCATABLE :: uP(:)
      CONTAINS
!        Initializes by partitioning the mesh and faces
         PROCEDURE :: ini => iniMsh 
!        Reads a mesh in parallel
         PROCEDURE :: read => readMsh
!        To deallocate the mesh
         PROCEDURE :: free => freeMsh
!        Partitioning a face
         PROCEDURE :: partFa => partFaMsh
      END TYPE mshType

      INTERFACE mshType
         PROCEDURE :: newMsh, newMshFl
      END INTERFACE mshType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
 
      FUNCTION newMsh(fName, nFa) RESULT(msh)
      TYPE(mshType) :: msh
      CHARACTER(LEN=*), INTENT(IN) :: fName
      INTEGER, INTENT(IN) :: nFa
      
      msh%name = ADJUSTL(fName)
      msh%nFa  = nFa
      ALLOCATE(msh%fa(nFa), msh%ub(nsd), msh%lb(nsd))

      RETURN
      END FUNCTION newMsh
!---------------------------------------------------------------------
      FUNCTION newMshFl(lst, fName) RESULT(msh)
      TYPE(lstType), INTENT(INOUT) :: lst
      CHARACTER(LEN=*), INTENT(IN) :: fName
      TYPE(mshType), TARGET :: msh

      INTEGER iFa, nFa
      REAL(KIND=8) scaleF
      CHARACTER(LEN=stdL) stmp,typ
      TYPE(fileType) f
      TYPE(lstType), POINTER :: lPtr, lPBC, lPty

      nFa  = lst%srch("Add face")
      io%o = " Number of available faces: "//nFa
      msh  = newMsh(fName,nFa)

      lPtr => lst%get(f,"vtk file path",1)
      CALL msh%read(file=f)
         
      scaleF = 1D0
      lPtr => lst%get(scaleF,"Mesh scale factor",lb=0D0)
      IF (ASSOCIATED(lPtr)) msh%x = msh%x*scaleF

      DO iFa=1, nFa
         lPBC => lst%get(stmp,"Add face",iFa)
         lPtr => lPBC%get(f,"vtk file path")
         lPty => lPBC%get(typ,"Face type")
         SELECT CASE (LOWER(TRIM(typ)))
            CASE ("inlet")
               msh%fa(iFa)%typ = 1
            CASE("outlet")
               msh%fa(iFa)%typ = 2
            CASE("wall")
               msh%fa(iFa)%typ = 3
            CASE DEFAULT
         io%e = "Select inlet, outlet, or wall Face type"         
         END SELECT
         CALL msh%readFace(iFa,stmp,file=f)
      END DO

      RETURN
      END FUNCTION newMshFl
!---------------------------------------------------------------------
      SUBROUTINE readMsh(lM, path, file)
      CLASS(mshType), INTENT(INOUT) :: lM
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: path
      TYPE(fileType), INTENT(IN), OPTIONAL :: file

      IF (cm%mas()) CALL readPMsh(lM, path, file)

!     Sending scalar data first
      CALL cm%bcast(lM%nEl)
      CALL cm%bcast(lM%nNo)
      CALL cm%bcast(lM%eType)
      CALL cm%bcast(lM%vol)
     
!     Setting global variables
      lM%gnEl = lM%nEl
      lM%gnNo = lM%nNo

      IF (cm%slv()) lM%eleType = eleType(lM%eType)

      RETURN
      END SUBROUTINE readMsh
!---------------------------------------------------------------------
!     This is for partitioning a single mesh
      SUBROUTINE iniMsh(lM, wgt, dof, guP, gmtl)
      CLASS(mshType), INTENT(INOUT) :: lM
      REAL(KIND=8), INTENT(IN), OPTIONAL :: wgt(:)
      INTEGER, INTENT(INOUT), OPTIONAL :: dof
      INTEGER, INTENT(INOUT), ALLOCATABLE, OPTIONAL :: guP(:)
      INTEGER, INTENT(INOUT), OPTIONAL :: gmtl(:)

      INTEGER, PARAMETER :: fid=1
      LOGICAL flag
      INTEGER p, a, Ac, e, Ec, edgecut, nNo, eNoNb
      CHARACTER(LEN=stdL) ctmp
      REAL, ALLOCATABLE :: fWgt(:)
      INTEGER, ALLOCATABLE :: part(:), gPart(:),tIEN(:,:), sCount(:), 
     2   disp(:)
      REAL(KIND=8), ALLOCATABLE :: tX(:,:)
     
      IF (lM%izd) THEN
         io%w = "Mesh <"//TRIM(lM%name)//"> already initialzied"
         RETURN
      END IF
      lM%izd = .TRUE.

      IF (cm%seq()) THEN
         ALLOCATE(lM%eD(0:1), lM%gN(lM%nNo), lM%gtl(lM%nNo),
     2      lM%otn(lM%nEl))
         ALLOCATE(lM%gIEN,SOURCE=lM%IEN)
         lM%eD(0) = 0
         lM%eD(1) = lM%nEl
         lM%gN    = (/(a, a=1, lM%nNo)/)
         lM%gtl   = (/(a, a=1, lM%nNo)/)
         lM%otn   = (/(a, a=1, lM%nEl)/)
         RETURN
      END IF
 
      IF (cm%slv()) THEN
         ALLOCATE(lM%gIEN(0,0))
      ELSE
         CALL MOVE_ALLOC(lM%IEN,lM%gIEN)
      END IF

      ALLOCATE(fwgt(cm%np()))
      fwgt = 1.0/REAL(cm%np())
      IF (PRESENT(wgt)) fwgt = REAL(wgt)
      
      ALLOCATE(lM%eD(0:cm%np()))
!     A draft of splitting the mesh between processors
!     lM%eD(i) represents first element which belong to cm%id()=i
      DO p=0, cm%np()
         lM%eD(p) = NINT(SUM(fwgt(1:p))*lM%gnEl)
      END DO
      lM%eD(cm%np()) = lM%gnEl

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      DO p=1, cm%np()
         disp(p)   = lM%eD(p-1)*lM%eNoN
         sCount(p) = lM%eD(p)*lM%eNoN - disp(p)
      END DO

      lM%nEl = lM%eD(cm%id() + 1) - lM%eD(cm%id())
      ALLOCATE(part(lM%nEl))

      ctmp = TRIM(appPath)//".part_"//TRIM(lM%name)//"_"//STR(cm%np())
     2   //".bin"
      INQUIRE(FILE=ctmp, EXIST=flag)
      IF (lM%eType.EQ.eType_NRB .OR. ANY(fwgt.EQ.1.0)) THEN
         part = cm%id()
      ELSE IF (flag) THEN
         OPEN(fid,FILE=ctmp, ACCESS='STREAM')
         CALL cm%read(part,fid)
         CLOSE(fid)
      ELSE 
         ALLOCATE(lM%IEN(lM%eNoN,lM%nEl))
!     Scattering the lM%gIEN array to processors
         CALL cm%scatter(lM%gIEN,lM%IEN)

!     This is to get eNoNb   
         SELECT CASE (lM%eType)
         CASE(eType_BRK) 
            eNoNb = 4
         CASE(eType_TET) 
            eNoNb = 3
         CASE(eType_WDG) 
            eNoNb = 3
         CASE(eType_TRI) 
            eNoNb = 2
         CASE(eType_BIL) 
            eNoNb = 2
         CASE(eType_BIQ) 
            eNoNb = 3
         CASE DEFAULT 
            io%e = "Undefined element type"
         END SELECT

!     The output of this call is "part" array in which part(i) is the
!     processor ID that element "i" belongs to
!     Doing partitioning, using ParMetis
         CALL SPLIT(lM%nEl, lM%eNoN, eNoNb, lM%IEN, cm%np(), lM%eD,
     2      fwgt, part, edgecut)
         IF (edgecut .EQ. 0) THEN
            io%w = " ParMETIS failed to partition the mesh"
            part = cm%id()
         ELSE IF (edgecut .GT. 0) THEN
            io%o = " ParMETIS partitioned the mesh by cutting "//
     2         STR(edgecut)//" elements"
!     LT 0 is for the case that all elements reside in one processor
         END IF
         DEALLOCATE(lM%IEN)
         OPEN(fid,FILE=ctmp, ACCESS='STREAM')
         CALL cm%write(part,fid)
         CLOSE(fid)
      END IF

!     Gathering the parts inside master.
!     gpart is a global version of part in which processor p = gpart(e)
!     is the owner of element "e"
      gPart = cm%gather(part)

      DEALLOCATE(part)
      ALLOCATE(lM%otn(lM%gnEl))
      IF (cm%mas()) THEN
         sCount = 0
         DO e=1, lM%gnEl
            p = gPart(e) + 1
            sCount(p) = sCount(p) + 1
         END DO
         DO p=1, cm%np()
            lM%eD(p) = lM%eD(p-1) + sCount(p)
         END DO

         ALLOCATE(tIEN(lM%eNoN,lM%gnEl))
!     Making the lM%IEN array in order, based on the cm%id() number in the 
!     master
         DO e=1, lM%gnEl
            p          = gPart(e)
            lM%eD(p)   = lM%eD(p) + 1
            Ec         = lM%eD(p)
            tIEN(:,Ec) = lM%gIEN(:,e)
            lM%otn(e)  = Ec
         END DO
         lM%gIEN = tIEN
         lM%eD(0) = 0
         DO p=1, cm%np()
            lM%eD(p) = lM%eD(p-1) + sCount(p)
         END DO
      END IF
      DEALLOCATE(gPart)
         
      CALL cm%bcast(lM%otn)
      CALL cm%bcast(lM%eD)
      
      lM%nEl = lM%eD(cm%id() + 1) - lM%eD(cm%id())
      ALLOCATE(lM%IEN(lM%eNoN,lM%nEl))

!     Now scattering the sorted lM%IEN to all processors
      IF (.NOT.ALLOCATED(tIEN)) ALLOCATE(tIEN(0,0))
      CALL cm%scatter(tIEN, lM%IEN)
      DEALLOCATE(tIEN)
     
!     Constructing the initial global to local pointer
!     lM%IEN: eNoN,nEl --> gnNo
!     lM%gtl: gnNo     --> nNo
!     lM%IEN: eNoN,nEl --> nNo
!     lM%gN:  nNo      --> gnNo
      ALLOCATE(lM%gtl(lM%gnNo))
      nNo    = 0
      lM%gtl = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            IF (lM%gtl(Ac) .EQ. 0) THEN
               nNo = nNo + 1
               lM%gtl(Ac) = nNo
            END IF
            lM%IEN(a,e) = lM%gtl(Ac)
         END DO
      END DO
      lM%nNo = nNo
      ALLOCATE(lM%gN(lM%nNo), tX(nsd,lM%nNo))
      DO Ac=1, lM%gnNo
         IF (lM%gtl(Ac) .NE. 0) lM%gN(lM%gtl(Ac)) = Ac
      END DO

!     Distributing X to processors
      CALL cm%local(tX, lM%x, lM%gN)
      IF (ALLOCATED(lM%x)) DEALLOCATE(lM%x)
      CALL MOVE_ALLOC(tX,lM%x)

!     Resetting lb and ub based on the size of the partitioned mesh
      lM%ub = MAXVAL(lM%x,2)
      lM%lb = MINVAL(lM%x,2)

!     Assigning unknown pointers if gmtl is provided
!     If things are not that simple mapping and converting other 
!     parameters. I will use an upper bound for guP, since there 
!     can be repeated nodes. aU is a container for QSORT.
!     part:  nNo  --> gdof
!     guP:   dof  --> gdof   (gPart is its tmp array)
!     uP:    gnNo --> gdof   (master-IN)
!     uP:    nNo  --> dof    (all-OUT)
!     gmtl:  gdof --> dof
      IF (PRESENT(gmtl)) THEN
         IF (.NOT.PRESENT(guP)) io%e = "msh%iniMsh: guP must be "//
     2      "present if gmtl is present"
         nNo  = lM%nNo
         ALLOCATE(part(nNo))
         CALL cm%local(part, lM%uP, lM%gN)
         IF (cm%mas()) DEALLOCATE(lM%uP)
         ALLOCATE(gPart(dof+nNo), lM%uP(nNo))
         DO a=1, dof
            Ac       = guP(a)
            gPart(a) = Ac
            gmtl(Ac) = a
         END DO
         DO a=1, nNo
            Ac = part(a)
            IF (gmtl(Ac) .EQ. 0) THEN
               dof        = dof + 1
               gmtl(Ac)   = dof
               lM%uP(a)   = dof
               gPart(dof) = Ac
            ELSE
               lM%uP(a) = gmtl(Ac)
            END IF
         END DO
         IF (ALLOCATED(guP)) DEALLOCATE(guP)
         ALLOCATE(guP(dof))
         guP = gPart(1:dof)
         DEALLOCATE(part, gPart)
      END IF

      RETURN
      END SUBROUTINE iniMsh
!---------------------------------------------------------------------
!     This routine partitions the face based on the already 
!     partitioned mesh
      SUBROUTINE partFaMsh(lM, iFa, gmtl)
      CLASS(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: iFa, gmtl(:)

      INTEGER e, a, Ac, Ec, i, j, nRing
      TYPE(faceType) :: gFa
      ASSOCIATE(fa => lM%fa(iFa))

!     A copy of the original face 
      gFa = fa
      DEALLOCATE(fa%gE, fa%IEN, fa%gN, fa%nV, fa%ring)
      IF (cm%mas()) DEALLOCATE(fa%uP)
      IF (cm%slv()) ALLOCATE(gFa%uP(gFa%nNo))
      CALL cm%bcast(gFa%uP)

!     Finding the number of lM%fas to allocate required space, also
!     maping global element number to processor element number
      fa%nEl = 0
      DO e=1, gFa%nEl
         Ec = lM%otn(gFa%gE(e))
         gFa%gE(e) = Ec
         IF (Ec.LE.lM%eD(cm%id()+1) .AND. Ec.GT.lM%eD(cm%id())) 
     2      fa%nEl = fa%nEl + 1
      END DO
      ALLOCATE(fa%gE(fa%nEl), fa%IEN(fa%eNoN,fa%nEl))
      fa%nNo = 0
      fa%dof = 0
      nRing  = 0 
      DO a=1, gFa%nNo
         Ac = gmtl(gFa%uP(a))
         IF (Ac .NE. 0) fa%dof = fa%dof + 1
         Ac = lM%gtl(gFa%gN(a))
         IF (Ac .NE. 0) fa%nNo = fa%nNo + 1
      END DO
      DO a=1, SIZE(gFa%ring)
         Ac = lM%gtl(gFa%ring(a))
         IF (Ac .NE. 0) nRing  = nRing + 1
      END DO
      ALLOCATE(fa%gN(fa%nNo), fa%uP(fa%dof), fa%nV(nsd,fa%nNo),
     2   fa%ring(nRing))
 
!     Time to form "face" structure in each processor
!     Copying the nodes which belong to this processor
      i = 0
      j = 0
      DO a=1, gFa%nNo
         Ac = gmtl(gFa%uP(a))
         IF (Ac .NE. 0) THEN
            j        = j + 1
            fa%uP(j) = Ac
         END IF
         Ac = lM%gtl(gFa%gN(a))
         IF (Ac .NE. 0) THEN
            i          = i + 1
            fa%gN(i)   = Ac
            fa%nV(:,i) = gFa%nV(:,a)
         END IF
      END DO
      i = 0
      DO a=1, SIZE(gFa%ring)
         Ac = lM%gtl(gFa%ring(a))
         IF (Ac .NE. 0) THEN
            i          = i + 1
            fa%ring(i) = Ac
         END IF
      END DO

!     And copying the element which belong to this processors
      j = 0
      DO e=1, gFa%nEl
         Ec = gFa%gE(e)
         IF (Ec.LE.lM%eD(cm%id()+1).AND.Ec.GT.lM%eD(cm%id())) THEN
            j = j + 1
            fa%gE(j) = Ec - lM%eD(cm%id())
            DO a=1, fa%eNoN
               fa%IEN(a,j) = lM%gtl(gFa%IEN(a,e))
            END DO
         END IF
      END DO
      CALL gFa%free()
      
      RETURN
      END ASSOCIATE
      END SUBROUTINE partFaMsh
!---------------------------------------------------------------------
      SUBROUTINE freeMsh(this)
      CLASS(mshType) :: this

      this%izd = .FALSE.
      IF (ALLOCATED(this%eD))   DEALLOCATE(this%eD)
      IF (ALLOCATED(this%gIEN)) DEALLOCATE(this%gIEN)
      IF (ALLOCATED(this%gN))   DEALLOCATE(this%gN)
      IF (ALLOCATED(this%uP))   DEALLOCATE(this%uP)
      CALL freePMsh(this)
      
      RETURN 
      END SUBROUTINE freeMsh

!#####################################################################
 
      SUBROUTINE CREATE_DUMMY_MESH(fName)
      CHARACTER(LEN=*), INTENT(IN) :: fName

      INTEGER, PARAMETER :: fid = 987
      INTEGER i

      OPEN(fid,FILE=TRIM(fName)//'.vtk')
      WRITE(fid,'(A)') "# vtk DataFile Version 2.0"
      WRITE(fid,'(A)') "Unstructured Grid"
      WRITE(fid,'(A)') "ASCII"
      WRITE(fid,'(A)') "DATASET UNSTRUCTURED_GRID"
      WRITE(fid,'(A)') "POINTS 8 double"
      WRITE(fid,'(A)') "1 1 1"
      WRITE(fid,'(A)') "0 1 0"
      WRITE(fid,'(A)') "1 0 0"
      WRITE(fid,'(A)') "0 1 1"
      WRITE(fid,'(A)') "0 0 1"
      WRITE(fid,'(A)') "1 0 1"
      WRITE(fid,'(A)') "1 1 0"
      WRITE(fid,'(A)') "0 0 0"
      WRITE(fid,'(A)') "CELLS 6 30"
      WRITE(fid,'(A)') "4     1    3    4    0"
      WRITE(fid,'(A)') "4     7    1    4    0"
      WRITE(fid,'(A)') "4     7    1    0    2"
      WRITE(fid,'(A)') "4     0    1    6    2"
      WRITE(fid,'(A)') "4     4    7    0    2"
      WRITE(fid,'(A)') "4     5    4    0    2"
      WRITE(fid,'(A)') "CELL_TYPES 6"
      DO i=1, 6 
         WRITE(fid,'(A)') "10"
      END DO
      CLOSE(fid)
      OPEN(fid,FILE=TRIM(fName)//'pFace.vtk')
      WRITE(fid,'(A)') "# vtk DataFile Version 3.0"
      WRITE(fid,'(A)') "vtk output"
      WRITE(fid,'(A)') "ASCII"
      WRITE(fid,'(A)') "DATASET POLYDATA"
      WRITE(fid,'(A)') "POINTS 4 double"
      WRITE(fid,'(A)') "1 1 1"
      WRITE(fid,'(A)') "0 1 0"
      WRITE(fid,'(A)') "0 1 1"
      WRITE(fid,'(A)') "1 1 0"
      WRITE(fid,'(A)') "POLYGONS 2           8"
      WRITE(fid,'(A)') "3       1       3       0"
      WRITE(fid,'(A)') "3       1       0       2"
      WRITE(fid,'(A)') "POINT_DATA 4"
      WRITE(fid,'(A)') "SCALARS scalars float"
      WRITE(fid,'(A)') "LOOKUP_TABLE default"
      WRITE(fid,'(A)') "1 2 4 7"
      CLOSE(fid)
      OPEN(fid,FILE=TRIM(fName)//'Face.vtk')
      WRITE(fid,'(A)') "# vtk DataFile Version 3.0"
      WRITE(fid,'(A)') "vtk output"
      WRITE(fid,'(A)') "ASCII"
      WRITE(fid,'(A)') "DATASET POLYDATA"
      WRITE(fid,'(A)') "POINTS 4 double"
      WRITE(fid,'(A)') "1 0 1"
      WRITE(fid,'(A)') "0 0 0"
      WRITE(fid,'(A)') "0 0 1"
      WRITE(fid,'(A)') "1 0 0"
      WRITE(fid,'(A)') "POLYGONS 2           8"
      WRITE(fid,'(A)') "3       2       1       3"
      WRITE(fid,'(A)') "3       0       2       3"
      WRITE(fid,'(A)') "POINT_DATA 4"
      WRITE(fid,'(A)') "SCALARS scalars float"
      WRITE(fid,'(A)') "LOOKUP_TABLE default"
      WRITE(fid,'(A)') "6 8 5 3"
      CLOSE(fid)
      OPEN(fid,FILE=TRIM(fName)//'nFace.vtk')
      WRITE(fid,'(A)') "# vtk DataFile Version 3.0"
      WRITE(fid,'(A)') "vtk output"
      WRITE(fid,'(A)') "ASCII"
      WRITE(fid,'(A)') "DATASET POLYDATA"
      WRITE(fid,'(A)') "POINTS 4 double"
      WRITE(fid,'(A)') "1 1 0"
      WRITE(fid,'(A)') "0 1 0"
      WRITE(fid,'(A)') "0 0 0"
      WRITE(fid,'(A)') "1 0 0"
      WRITE(fid,'(A)') "POLYGONS 2           8"
      WRITE(fid,'(A)') "3       2       1       3"
      WRITE(fid,'(A)') "3       0       1       3"
      WRITE(fid,'(A)') "POINT_DATA 4"
      WRITE(fid,'(A)') "SCALARS scalars float"
      WRITE(fid,'(A)') "LOOKUP_TABLE default"
      WRITE(fid,'(A)') "7 2 8 3"
      CLOSE(fid)

      RETURN
      END SUBROUTINE CREATE_DUMMY_MESH
!---------------------------------------------------------------------
!     To test the current implementation.
      SUBROUTINE TEST_MSHMOD()

      TYPE(mshType) :: msh
      REAL(KIND=8), ALLOCATABLE :: gX(:,:)
 
!     Writing a dummy mesh to perform testing
      IF (cm%mas()) CALL CREATE_DUMMY_MESH('.dummy')
      msh = mshType("My_mesh",2) !mesh_name, num_of_face,comu
      IF (msh%nFa.NE.2) io%e = "Issue with newMsh"
      io%o = "newMsh: "//CLR("(PASSED)",3)
      CALL msh%read(".dummy.vtk") !Reading mesh from vtk file
      io%o = "msh%read: "//CLR("(PASSED)",3)
      IF (msh%nEl.NE.6 .OR. msh%nNo.NE.8 .OR. .NOT.ISZERO(msh%vol,1D0)) 
     2   io%e = "Issue with pMsh%n*"
      io%o = "msh%n*: "//CLR("(PASSED)",3)
      IF (cm%mas()) THEN
         IF (ANY(msh%IEN(:,1).NE.(/4,2,5,1/))) io%e = "Issue with IEN"
         io%o = "pMsh%IEN: "//CLR("(PASSED)",3)
         IF (ANY(msh%x(:,1).NE.1D0)) io%e = "Issue with x"
         io%o = "pMsh%x: "//CLR("(PASSED)",3)
      END IF
      CALL msh%readFace(1,"Per_face",".dummypFace.vtk")! ID, name, path
      CALL msh%readFace(2,"My_face",".dummyFace.vtk") ! ID, name, path
      IF (ANY(msh%fa(1)%IEN(:,1).NE.(/2,7,1/)) .OR. 
     2   msh%fa(2)%name.NE."My_face") io%e = "Issue with face%IEN"
      io%o = "msh%readFace: "//CLR("(PASSED)",3)
      ALLOCATE(gX,SOURCE=msh%x)
      IF (.NOT.ALLOCATED(gX)) ALLOCATE(gX(0,0))
      CALL msh%ini() !Partitioning the mesh
      IF (cm%np().EQ.2 .AND. ANY(msh%eD.NE.(/0,3,6/))) 
     2   io%e = "Issue with ini"
      IF (cm%seq() .AND. ANY(msh%eD.NE.(/0,6/))) 
     2   io%e = "Issue with ini"
      io%o = "msh%ini: "//CLR("(PASSED)",3)
      IF (.NOT.ISZERO(msh%fa(1)%area,1d0)) io%e = "Issue with setupFace"
      io%o = "face%setup: "//CLR("(PASSED)",3)
      
      RETURN
      END SUBROUTINE TEST_MSHMOD

      END MODULE MSHMOD

