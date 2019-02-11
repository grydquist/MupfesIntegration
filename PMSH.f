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
      MODULE PMSHMOD
      USE ELEMOD
      USE LSTMOD
      IMPLICIT NONE

!     The face type containing mesh at boundary
      TYPE, EXTENDS(eleType) :: faceType
c         PRIVATE
!        Number of elements
         INTEGER :: nEl = 0
!        Number of nodes
         INTEGER :: nNo = 0
!        Number if degrees of freedom on this face
         INTEGER :: dof = 0
!        Global element Ids [fa%nEl --> msh%nEl]
         INTEGER, ALLOCATABLE :: gE(:)
!        Global node Ids [fa%nNo --> msh%nNo]
         INTEGER, ALLOCATABLE :: gN(:)
!        Pointer to unknown ID [fa%dof --> dmn%dof]
         INTEGER, ALLOCATABLE :: uP(:)
!        List of nodes on the perimeter of face (master)
         INTEGER, ALLOCATABLE :: ring(:)
!        Connectivity array [fa%eNoN,fa%nEl --> msh%nNo]
         INTEGER, ALLOCATABLE :: IEN(:,:)
!        Surface area
         REAL(KIND=8) area
!        Normal vector to each nodal point [nsd x nNo]
         REAL(KIND=8), ALLOCATABLE :: nV(:,:)
!        Face name for flux files
         CHARACTER(LEN=stdL) :: name = DEFAULT_NAME
!        Type of face (1 = inlet, 2 = outlet, 3 = wall)
         INTEGER :: typ
      CONTAINS
!        To deallocate the face
         PROCEDURE :: free => freeFace
      END TYPE faceType
!--------------------------------------------------------------------
!     Parent mesh type (single processor)
      TYPE, EXTENDS(eleType), ABSTRACT :: pMshType
c         PRIVATE
!        Number of elements (knot spanes)
         INTEGER :: nEl = 0
!        Number of nodes
         INTEGER :: nNo = 0
!        Number of faces         
         INTEGER :: nFa = 0
!        Volume
         REAL(KIND=8) vol
!        The connectivity array [eNoN,nEl --> nNo]
         INTEGER, ALLOCATABLE :: IEN(:,:)
!        Lower bound of the mesh [nsd]
         REAL(KIND=8), ALLOCATABLE :: lb(:)
!        Upper bound of the mesh [nsd]
         REAL(KIND=8), ALLOCATABLE :: ub(:)
!        Position of nodes [nsd x nNo]
         REAL(KIND=8), ALLOCATABLE :: x(:,:)
!        Mesh name
         CHARACTER(LEN=stdL) :: name = DEFAULT_NAME
!        Faces are stored in this variable [nFa]
         TYPE(faceType), ALLOCATABLE :: fa(:)
      CONTAINS
!        To integrate a variable over the volume
         PROCEDURE :: integSPMsh
         PROCEDURE :: integGPMsh
         PROCEDURE :: integSFacePMsh
         PROCEDURE :: integVFacePMsh
         PROCEDURE :: integGFacePMsh
         GENERIC :: integ => integSFacePMsh, integVFacePMsh, 
     2      integGFacePMsh, integSPMsh, integGPMsh
!        To deallocate the mesh
         PROCEDURE :: free => freePMsh
!        Computes normal to a face at a given element and Gauss point
         PROCEDURE :: normal => normalPMsh
!        Finds nodes on two faces sharing the same position
         PROCEDURE :: matchFace => matchFacePMsh
!        Reads the face from a vtk file (to be combined with readMsh)
         PROCEDURE :: readFace => readFacePMsh
!        Finds the ring around a face
         PROCEDURE, PRIVATE :: ringFace => ringFacePMsh
!        Setup the face: computes normal and area
         PROCEDURE, PRIVATE :: setupFace => setupFacePMsh
!        Sets face EBC
         PROCEDURE, PRIVATE :: ebcFace => ebcFacePMsh
!        Sets face NBC
         PROCEDURE, PRIVATE :: nbcFace => nbcFacePMsh
!        Checks to make sure everything is in order
         PROCEDURE, PRIVATE :: check => checkPMsh
!        Reads the mesh from a vtk file
         PROCEDURE, PRIVATE :: read => readPMsh
      END TYPE pMshType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
!     Reading face from a VTK file
      SUBROUTINE readFacePMsh(lM, iFa, name, path, file)
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: iFa
      CHARACTER(LEN=*), INTENT(IN) :: name
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: path
      TYPE(fileType), INTENT(IN), OPTIONAL :: file

      INTEGER e, a, fid
      CHARACTER(LEN=stdL) rLine
      TYPE(fileType) f
      INTEGER, ALLOCATABLE :: nbc(:), gnID(:)
      ASSOCIATE(fa => lM%fa(iFa))

      fa%name = ADJUSTL(name)
      IF (PRESENT(path)) THEN
         f = fileType(path)
      ELSE
         IF (.NOT.PRESENT(file)) io%e = 
     2      "readFace: path or file must be specified"
         f = file
      END IF
      CALL f%open('r')
      fid = f%id()
      DO
         READ(fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(1:10) .EQ. 'POINT_DATA') EXIT
      END DO
      rLine = rLine(11:)
      READ (rLine,*) a
      ALLOCATE(nbc(a))
      DO
         READ (fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
!     Continuing until to get to a block of nodal data
         IF (rLine(1:7) .EQ. 'SCALARS') THEN
            READ (fid,'(a)',END=001) rLine
            EXIT
         ELSE 
            CYCLE
         END IF
      END DO
!     Reading global ID of the nodes
      READ(fid,*) nbc

      REWIND(fid)
      DO
         READ(fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(1:8) .EQ. 'POLYGONS') EXIT
      END DO
      rLine = rLine(9:)
      READ(rLine,*) fa%nEl
      READ(fid,'(a)',END=001) rLine
      fa%eNoN = CheckNoNumbers(rLine) - 1
!     Finding element type
      fa%eType = ELEMENT_TYPE(nsd-1, fa%eNoN)
      fa%eleType = eleType(fa%eType)
      ALLOCATE(fa%gE(fa%nEl), fa%IEN(fa%eNoN,fa%nEl), 
     2   gnID(fa%eNoN))
      READ (rLine,*) a, gnID
      IF (a .NE. fa%eNoN) io%e = "Inconsistent number of columns"
      DO a=1, fa%eNoN
         fa%IEN(a,1) = nbc(gnID(a)+1)
      END DO
!     And reading the ebc data
      DO e=2, fa%nEl
         READ(fid,*,END=001) a, gnID
         DO a=1, fa%eNoN
            fa%IEN(a,e) = nbc(gnID(a)+1)
         END DO
      END DO
      CALL f%close()

      IF (cm%mas()) CALL lM%setupFace(iFa)
      IF (cm%mas()) a = SIZE(fa%ring)
      CALL cm%bcast(a)
      CALL cm%bcast(fa%nNo)
      CALL cm%bcast(fa%area)
      CALL cm%bcast(fa%gE)
      IF (cm%slv()) ALLOCATE(fa%gN(fa%nNo), fa%nV(nsd,fa%nNo),
     2   fa%ring(a))
      CALL cm%bcast(fa%ring)
      CALL cm%bcast(fa%gN)
      CALL cm%bcast(fa%nV)

      RETURN
 001  io%e = "A block of data is missing from the BC vtk file"
      END ASSOCIATE
      END SUBROUTINE readFacePMsh
!---------------------------------------------------------------------
!     Setting up the faces (ebc/nbc/nV). To be called by master only.     
      SUBROUTINE setupFacePMsh(lM, iFa)
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: iFa

      LOGICAL flag
      INTEGER e, a, Ac, g, nNo
      REAL(KIND=8) tmp, n(nsd)
      REAL(KIND=8), ALLOCATABLE :: s(:), sV(:,:)
      ASSOCIATE(fa => lM%fa(iFa))

!     And calculating enc, nbc and ring
      CALL lM%ebcFace(iFa)
      CALL lM%nbcFace(iFa)
      CALL lM%ringFace(iFa)

      nNo = lM%nNo
!     Calculating the center of the face, diameter and its area
      ALLOCATE(s(nNo), sV(nsd,nNo), fa%nV(nsd,fa%nNo))
      s       = 1D0
      sV      = 0D0
      fa%area = lM%integ(iFa,s)
!     Making sure area is not zero, since it will cause issues later on
      IF (ISZERO(fa%area)) io%w = "<"//TRIM(fa%name)//"> area is zero"
      DO e=1, fa%nEl
c         IF (fa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), fa, e)
         DO g=1, fa%nG
            n = lM%normal(iFa,e,g)
            DO a=1, fa%eNoN
               Ac       = fa%IEN(a,e)
               sV(:,Ac) = sV(:,Ac) + n*fa%N(a,g)*fa%w(g)
            END DO
         END DO
      END DO
      
      flag = .TRUE.
      DO a=1, fa%nNo
         Ac  = fa%gN(a)
         tmp = NORM2(sV(:,Ac))
         IF (ISZERO(tmp)) THEN
            IF (flag) THEN
               io%w = "Skipping normal calculation of node "//a//
     2            " in face <"//TRIM(fa%name)//">"
               flag = .FALSE.
            END IF
            fa%nV(:,a) = 0D0
            fa%nV(1,a) = 1D0
            CYCLE
         END IF
         fa%nV(:,a) = sV(:,Ac)/tmp
      END DO

      RETURN
      END ASSOCIATE
      END SUBROUTINE setupFacePMsh
!---------------------------------------------------------------------
!     Construct %gE based on IEN and face%IEN
      SUBROUTINE ebcFacePMsh(lM, iFa)
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: iFa

      INTEGER e, a, b, Ac, ie, maxnEtN
      INTEGER, ALLOCATABLE :: nbc(:), ndEs(:,:), bcNdEs(:,:)
      ASSOCIATE(fa => lM%fa(iFa))

!     First constructing an array, ndEs, with size
!     nNo*(max no element touching a node) that row "a" contains 
!     elements that have "a" as a node. Here "nbc" array keeps track of
!     the length ndEs
      ALLOCATE(nbc(lM%nNo))
      maxnEtN = 3*lM%eNoN**3
 003  maxnEtN = maxnEtN + 5
      IF (ALLOCATED(ndEs)) DEALLOCATE (ndEs)
      ALLOCATE (ndEs(lM%nNo,maxnEtN))
      nbc = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            nbc(Ac) = nbc(Ac) + 1
            IF (nbc(Ac) .GT. maxnEtN) THEN
               io%w = "Increasing maxnEtN to "//maxnEtN
               GOTO 003
            END IF
            ndEs(Ac,nbc(Ac)) = e
         END DO
      END DO
!     Now I am going to construct another array that contains the
!     elements that touch boundary nodes in the first index and the
!     number of repetition of those elements in the second index
      ALLOCATE(bcNdEs(maxnEtN*fa%eNoN,2))
      DO e=1, fa%nEl
         bcNdEs = 0
         DO a=1, fa%eNoN
            Ac = fa%IEN(a,e)
            DO ie=1, nbc(Ac)
               DO b=1, maxnEtN*fa%eNoN
                  IF (bcNdEs(b,1) .EQ. 0) THEN
                     bcNdEs(b,1) = ndEs(Ac,ie)
                     bcNdEs(b,2) = 1
                     EXIT
                  ELSE IF (bcNdEs(b,1) .EQ. ndEs(Ac,ie)) THEN
                     bcNdEs(b,2) = bcNdEs(b,2) + 1
                     EXIT
                  END IF
               END DO
            END DO
         END DO
!     Now if the number of repetition of an element is equal to face 
!     eNoN, then that element should have a face shared with the 
!     boundary element
         DO b=1, maxnEtN*fa%eNoN
            IF (bcNdEs(b,2) .EQ. fa%eNoN) THEN
               fa%gE(e) = bcNdEs(b,1)
               EXIT
            ELSE IF (bcNdEs(b,2) .EQ. 0) THEN
               io%e = "global element number for element <"//e
     2            //"> of face <"//TRIM(fa%name)//"> was not found"
            END IF
         END DO
      END DO
      
      RETURN
      END ASSOCIATE
      END SUBROUTINE ebcFacePMsh
!---------------------------------------------------------------------
!     Making sure all the faces contain in range values
      SUBROUTINE nbcFacePMsh(lM, iFa)
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: iFa

      INTEGER Ac, e, Ec, a
      INTEGER, ALLOCATABLE :: incNd(:)
      ASSOCIATE(fa => lM%fa(iFa))

!     Finding the nodes that belongs to this face and counting the
!     number of the nodes
      ALLOCATE(incNd(lM%nNo))
      incNd   = 0
      fa%nNo = 0
      DO e=1, fa%nEl
         Ec = fa%gE(e)
         IF (Ec.GT.lM%nEl .OR. Ec.LE.0) THEN
            io%e = "Global ID of element "//e//" of face <"
     2         //TRIM(fa%name)//"> is out of range"
         END IF
         DO a=1, fa%eNoN
            Ac = fa%IEN(a,e)
            IF (Ac.GT.lM%nNo .OR. Ac.LE.0) THEN
               io%e = "IEN of element "//e//" of face <"
     2            //TRIM(fa%name)//"> is out of range"
            END IF
            IF (ALL(lM%IEN(:,Ec) .NE. Ac)) THEN
               io%e = "Parent element of "//e//" of face <"
     2            //TRIM(fa%name)//"> does not contain same"
     3            //" nodes as the boundary"
            END IF
            IF (incNd(Ac) .EQ. 0) THEN
               fa%nNo = fa%nNo + 1
               incNd(Ac) = 1
            END IF
         END DO
      END DO
        
      ALLOCATE(fa%gN(fa%nNo))
      a = 0
      DO Ac=1, lM%nNo
         IF (incNd(Ac) .NE. 0) THEN
            a = a + 1
            fa%gN(a) = Ac
         END IF
      END DO

      RETURN
      END ASSOCIATE
      END SUBROUTINE nbcFacePMsh
!---------------------------------------------------------------------
!     Finds the ring around a face
      SUBROUTINE ringFacePMsh(lM, iFa)
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: iFa

      INTEGER e, a, b, Ac, uel, i, j, ei, ej, c, nnof
      LOGICAL, ALLOCATABLE :: inc(:)
      INTEGER, ALLOCATABLE :: tmp(:,:), gtl(:), adj(:,:)
      ASSOCIATE(fa => lM%fa(iFa))

      ALLOCATE(gtl(lM%nNo), adj(fa%nEf,fa%nEl), inc(fa%nNo))
      DO a=1, fa%nNo
         gtl(fa%gN(a)) = a
      END DO
!     Number of nodes on the boundary of face element
      IF (MOD(fa%eNoN,fa%nEf) .NE. 0) io%e = "ringFace: unexpected"//
     2   " combination of fa%eNoN and fa%nEf"
      nnof = fa%eNoN/fa%nEf + 1

      uel = 100
 001  uel = NINT(1.5*REAL(uel))
      ALLOCATE(tmp(uel,fa%nNo))
      tmp = 0
      DO e=1, fa%nEl
!     Including this element to tmp list
         DO a=1, fa%eNoN
            Ac = gtl(fa%IEN(a,e))
            i  = 1
            DO WHILE(tmp(i,Ac).NE.0)
               i = i + 1
!      A zero at the end is needed, hence GE
               IF (i .GE. uel) THEN
                  DEALLOCATE(tmp)
                  io%w = "ringFace: increasing uel, mesh might be bad"
                  GOTO 001
               END IF
            END DO
            tmp(i,Ac) = e
         END DO
      END DO

      adj = 0
      DO Ac=1, fa%nNo
         i = 1
         DO WHILE (tmp(i,Ac) .NE. 0)
            ei = tmp(i,Ac)
            j  = i + 1
            DO WHILE (tmp(j,Ac) .NE. 0)
               ej = tmp(j,Ac)
!     We have elements i and j that at least share one node. Lets see if
!     they share a second node
               c = 0
               DO a=1, fa%eNoN
                  DO b=1, fa%eNoN
                     IF (fa%IEN(a,ei) .EQ. fa%IEN(b,ej)) c = c + 1
                  END DO
               END DO
               IF (c .GT. nnof) io%e = "ringFace: too many nnof"
               IF (c .EQ. nnof) THEN
                  a = 1
                  DO WHILE (adj(a,ei).NE.0 .AND. adj(a,ei).NE.ej)
                     a = a + 1
                     IF (a .GT. fa%nEf) io%e = "ringFace: too many nei."
                  END DO
                  adj(a,ei) = ej
                  
                  a = 1
                  DO WHILE (adj(a,ej).NE.0 .AND. adj(a,ej).NE.ei)
                     a = a + 1
                     IF (a .GT. fa%nEf) io%e = "ringFace: too many nei."
                  END DO
                  adj(a,ej) = ei
               END IF
!     Search failed going to the next 
               j = j + 1
            END DO
            i = i + 1
         END DO
      END DO
      inc = .FALSE.
      DO e=1, fa%nEl
         IF (adj(fa%nEf,e) .NE. 0) CYCLE
!     Only considering nodes that are on the edge
         DO a=1, fa%nEf
            Ac = gtl(fa%IEN(a,e))
            c = 0
            i = 1
            DO WHILE (adj(i,e) .NE. 0)
               ei = adj(i,e)
               DO b=1, fa%nEf
                  IF (gtl(fa%IEN(b,ei)) .EQ. Ac) c = c + 1
               END DO
               i = i + 1
            END DO
            IF (c .LE. 1) inc(Ac) = .TRUE.
         END DO
      END DO

      c = COUNT(inc)
      ALLOCATE(fa%ring(c))
      c = 0
      DO Ac=1, fa%nNo
         IF (inc(Ac)) THEN
            c = c + 1
            fa%ring(c) = fa%gN(Ac)
         END IF
      END DO
      DEALLOCATE(tmp, inc, gtl, adj)

      RETURN
      END ASSOCIATE
      END SUBROUTINE ringFacePMsh
!---------------------------------------------------------------------
!     This routine returns a vector at element "e" and Gauss point 
!     "g" of face "fa" that is the normal weigthed by Jac, i.e. 
!     Jac = NORM2(n). 
      FUNCTION normalPMsh(lM, iFa, e, g) RESULT(nV)
      CLASS(pMshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: iFa, e, g
      REAL(KIND=8) :: nV(nsd)

      INTEGER a, Ac, i
      REAL(KIND=8) xXi(nsd,nsd-1), v(nsd), lX(nsd,lM%fa(iFa)%eNoN)
      ASSOCIATE(fa => lM%fa(iFa))

!     Correcting the position vector if mesh is moving
!     Also calculating surface deflation
      xXi = 0D0
      DO a=1, fa%eNoN
         Ac = fa%IEN(a,e)
         lX(:,a) = lM%x(:,Ac)
!         IF (mvMsh) lX(:,a) = lX(:,a) + Do(nsd+2:2*nsd+1,Ac)
         DO i=1, nsd-1
            xXi(:,i) = xXi(:,i) + fa%Nx(i,a,g)*lX(:,a)
         END DO
      END DO
      nV = CROSS(xXi)

!     Changing the sign if neccessary. First finding a node
!     outside of the face, in the parent element
      Ac = 1 ! To avoid warnings
      DO a=1, lM%eNoN
         Ac = lM%IEN(a,fa%gE(e))
         IF (ALL(fa%IEN(:,e).NE.Ac)) EXIT
      END DO
      IF (a .GT. lM%eNoN) io%e = 
     2   "normalPMsh: Unable to find a node outside of the face"

      v = lX(:,1) - lM%x(:,Ac)
      IF (DOT_PRODUCT(nV,v) .LT. 0D0) nV = -nV

      RETURN
      END ASSOCIATE
      END FUNCTION normalPMsh
!---------------------------------------------------------------------
!     This is match isoparameteric faces to each other. Project nodes 
!     from two adjacent meshes to each other based on a L2 norm.
      FUNCTION matchFacePMsh(iM, jM, iFa, iPFa) RESULT(lPrj)
      CLASS(pMshType), INTENT(IN) :: iM, jM
      INTEGER, INTENT(IN) :: iFa, iPFa
      TYPE(stackType) :: lPrj
 
      TYPE blkType
         INTEGER :: n = 0
         INTEGER, ALLOCATABLE :: gN(:)
      END TYPE

      LOGICAL nFlt(nsd)
      INTEGER nBkd, i, a, b, Ac, Bc, iBk, nBk
      REAL(KIND=8) v(nsd), xMin(nsd), xMax(nsd), dx(nsd)
      INTEGER, ALLOCATABLE :: nodeBlk(:)
      TYPE(blkType), ALLOCATABLE :: blk(:)
      ASSOCIATE(fa => iM%fa(iFa), pFa => jM%fa(iPFa))

      IF (fa%eType .NE. pFa%eType) 
     2   io%e = "Incompatible faces in MATCHFACES"

!     We want to have approximately 1000 nodes in each block. So we
!     calculate nBkd, which is the number of separate blockes in eahc
!     direction, based on that.
      a    = pFa%nNo
      nBkd = NINT((REAL(a,8)/1D3)**(3.33D-1))
      IF (nBkd .EQ. 0) nBkd = 1
      nBk  = nBkd**nsd
      ALLOCATE(nodeBlk(a), blk(nBk))

!     Finding the extends of the domain and size of each block
      DO i=1, nsd
         xMin(i) = MIN(MINVAL(iM%x(i,fa%gN)), MINVAL(jM%x(i,pFa%gN)))
         xMax(i) = MAX(MAXVAL(iM%x(i,fa%gN)), MAXVAL(jM%x(i,pFa%gN)))
         IF (xMin(i) .LT. 0D0) THEN
            xMin(i) = xMin(i)*(1D0+eps)
         ELSE
            xMin(i) = xMin(i)*(1D0-eps)
         END IF
         IF (xMax(i) .LT. 0D0) THEN
            xMax(i) = xMax(i)*(1D0-eps)
         ELSE
            xMax(i) = xMax(i)*(1D0+eps)
         END IF
      END DO
      dx = (xMax - xMin)/REAL(nBkd,8)
      nFlt = .TRUE.
      DO i=1, nsd
         IF (ISZERO(dx(i))) nFlt(i) = .FALSE.
      END DO

!     First finding an estimation for size of each block
      blk%n = 0
      DO a=1, pFa%nNo
         Ac  = pFa%gN(a)
         iBk = FINDBLK(jM%x(:,Ac))
         nodeBlk(a) = iBk
         blk(iBk)%n = blk(iBk)%n + 1
      END DO
      DO iBk=1, nBk
         ALLOCATE(blk(iBk)%gN(blk(iBk)%n))
      END DO
      blk%n = 0
      DO a=1, pFa%nNo
         Ac  = pFa%gN(a)
         iBk = nodeBlk(a)
         blk(iBk)%n = blk(iBk)%n + 1
         blk(iBk)%gN(blk(iBk)%n) = Ac
      END DO

!     Doing the calculation for every single node on this face
      DO a=1, fa%nNo
         Ac  = fa%gN(a)
         iBk = FINDBLK(iM%x(:,Ac))
!     Checking every single node on the other face
         DO b=1, blk(iBk)%n
            Bc = blk(iBk)%gN(b)
            IF (iM%name.EQ.jM%name .AND. Ac.EQ.Bc) CYCLE
            v  = jM%x(:,BC) - iM%x(:,Ac)
            IF (ISZERO(NORM2(v))) CALL lPrj%push((/Ac,Bc/))
         END DO
      END DO
      IF (.NOT.lPrj%pull(a)) io%e = "Unable to find any matching node"//
     2   " between <"//TRIM(fa%name)//"> and <"//TRIM(pFa%name)//">"
      CALL lPrj%push(a)

      RETURN
      END ASSOCIATE
      CONTAINS 
!`````````````````````````````````````````````````````````````````````
      INTEGER FUNCTION FINDBLK(x)
      REAL(KIND=8), INTENT(IN) :: x(nsd)

      INTEGER i, j, k

      i = 1 
      j = 1
      k = 1
      IF (nFlt(1)) i = INT((x(1) - xMin(1))/dx(1)) 
      IF (nFlt(2)) j = INT((x(2) - xMin(2))/dx(2)) 
      IF (i .EQ. nBkd) i = nBkd - 1
      IF (j .EQ. nBkd) j = nBkd - 1
      IF (nsd .EQ. 3) THEN
         IF (nFlt(3)) k = INT((x(3) - xMin(3))/dx(3)) 
         IF (k .EQ. nBkd) k = nBkd - 1
         FINDBLK = k + (j + i*nBkd)*nBkd + 1
      ELSE ! nsd .EQ. 2
         FINDBLK = j + i*nBkd + 1
      END IF

      RETURN
      END FUNCTION FINDBLK
      END FUNCTION matchFacePMsh
!---------------------------------------------------------------------
!     These set of routines destroy an object.
      SUBROUTINE freeFace(this)
      CLASS(faceType) :: this

      this%nEl = 0
      this%nNo = 0
      this%dof = 0

      IF (ALLOCATED(this%IEN))  DEALLOCATE(this%IEN)
      IF (ALLOCATED(this%nV))   DEALLOCATE(this%nV)
      IF (ALLOCATED(this%gN))   DEALLOCATE(this%gN)
      IF (ALLOCATED(this%gE))   DEALLOCATE(this%gE)
      IF (ALLOCATED(this%uP))   DEALLOCATE(this%uP)
      IF (ALLOCATED(this%ring)) DEALLOCATE(this%ring)

      RETURN 
      END SUBROUTINE freeFace
!---------------------------------------------------------------------
!     This routine integrate u over the surface faId.
      FUNCTION integSFacePMsh(lM, iFa, u) RESULT(s)
      CLASS(pMshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: iFa
      REAL(KIND=8), INTENT(IN) :: u(:)
      REAL(KIND=8) s

      INTEGER a, e, g, Ac
      REAL(KIND=8) uHat, Jac, n(nsd)
      ASSOCIATE(fa => lM%fa(iFa))

      IF (SIZE(u) .NE. lM%nNo) io%e = "integSFace: SIZE(u).NE.lM%nNo"
      s = 0D0
      DO e=1, fa%nEl
!     Updating the shape functions, if this is a NURB      
c         IF (fa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(fa%iM), fa, e)
         DO g=1, fa%nG
            n   = lM%normal(iFa, e, g)
            Jac = NORM2(n)
!     Calculating the function value
            uHat = 0D0
            DO a=1, fa%eNoN
               Ac   = fa%IEN(a,e)
               uHat = uHat + u(Ac)*fa%N(a,g)
            END DO
!     Now integrating
            s = s + Jac*fa%w(g)*uHat
         END DO
      END DO

      RETURN
      END ASSOCIATE
      END FUNCTION integSFacePMsh
!---------------------------------------------------------------------
!     This routine integrate u over the surface faId.
      FUNCTION integVFacePMsh(lM, iFa, u) RESULT(s)
      CLASS(pMshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: ifa
      REAL(KIND=8), INTENT(IN) :: u(:,:)
      REAL(KIND=8) s

      INTEGER a, i, e, Ac, g
      REAL(KIND=8) uHat, n(nsd)
      ASSOCIATE(fa => lM%fa(iFa))

      IF (SIZE(u,2).NE.lM%nNo) io%e = "integVFace: SIZE(u).NE.lM%nNo"
      s = 0D0
      DO e=1, fa%nEl
!     Updating the shape functions, if this is a NURB      
c         IF (fa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(fa%iM), fa, e)
         DO g=1, fa%nG
!     Since sHat must be multiplied by Jac, I don't normalize it here
            n   = lM%normal(iFa, e, g)
!     Calculating the function value
            uHat = 0D0
            DO a=1, fa%eNoN
               Ac = fa%IEN(a,e)
               DO i=1, nsd
                  uHat = uHat + fa%N(a,g)*u(i,Ac)*n(i)
               END DO
            END DO
!     Now integrating 
            s = s + fa%w(g)*uHat
         END DO
      END DO

      RETURN
      END ASSOCIATE
      END FUNCTION integVFacePMsh
!---------------------------------------------------------------------
!     This routine integrate s(lb:ub,:) over the surface faId.
      FUNCTION integGFacePMsh(lM, iFa, u, lb, ub) RESULT(s)
      CLASS(pMshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: iFa
      REAL(KIND=8), INTENT(IN) :: u(:,:)
      INTEGER, INTENT(IN) :: lb
      INTEGER, INTENT(IN), OPTIONAL :: ub
      REAL(KIND=8) s

      INTEGER a, fub
      REAL(KIND=8), ALLOCATABLE :: sclr(:), vec(:,:)
      ASSOCIATE(fa => lM%fa(iFa))

      fub = lb
      IF (PRESENT(ub)) fub = ub
      IF (SIZE(u,2).NE.lM%nNo) io%e = "integGFace: SIZE(u).NE.lM%nNo"

      s = 0D0
      IF (fub-lb+1 .EQ. nsd) THEN
         ALLOCATE(vec(nsd,lM%nNo))
         DO a=1, lM%nNo
            vec(:,a) = u(lb:fub,a)
         END DO
         s = lM%integVFacePMsh(iFa,vec)
      ELSE IF (lb .EQ. fub) THEN
         ALLOCATE (sclr(lM%nNo))
         DO a=1, lM%nNo
            sclr(a) = u(lb,a)
         END DO
         s = lM%integSFacePMsh(iFa,sclr)
      ELSE
         io%e = "IntegG: nexpected dof"
      END IF

      RETURN
      END ASSOCIATE
      END FUNCTION integGFacePMsh

!#####################################################################
!     This routine is desinged to read isoparametric meshes and 
!     retype it upon request
      SUBROUTINE readPMsh(lM, path, file)
      CLASS(pMshType), INTENT(INOUT) :: lM
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: path
      TYPE(fileType), INTENT(IN), OPTIONAL :: file

      INTEGER i, e, fid
      CHARACTER(LEN=stdL) stmp
      TYPE(fileType) f
      REAL(KIND=8), ALLOCATABLE :: rtmp(:), x(:,:)

      IF (PRESENT(path)) THEN
         f = fileType(path)
      ELSE
         IF (.NOT.PRESENT(file)) io%e = 
     2      "readMsh: path or file must be defined"
         f = file
      END IF
      CALL f%open('r')
      fid = f%id()
!     Reading vtk file type
      DO WHILE (.TRUE.)
         READ(fid,'(a)',END=001) stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(1:25) .EQ. 'DATASET UNSTRUCTURED_GRID') EXIT
      END DO
      READ(fid,'(a)',END=001) stmp
      stmp = ADJUSTL(stmp)
      stmp = stmp(7:)
      READ(stmp,*) lM%nNo
      ALLOCATE(lM%x(nsd,lM%nNo))
      IF (nsd .EQ. 2) THEN
         ALLOCATE(x(3,lM%nNo))
         READ(fid,*,END=001) x
         lM%x = x(1:2,:)
         DEALLOCATE(x)
      ELSE
         READ(fid,*,END=001) lM%x
      END IF

      DO WHILE (.TRUE.)
         READ(fid,'(a)',END=001) stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(1:5) .EQ. 'CELLS') EXIT
      END DO
      stmp = stmp(6:)
      READ (stmp,*) lM%nEl, i
      IF (MOD(i,lM%nEl) .NE. 0) io%e ="Faulty value in the vtk file"
      lM%eNoN  = i/lM%nEl - 1
      lM%eType = ELEMENT_TYPE(nsd, lM%eNoN)
      lM%eleType = eleType(lM%eType)
      ALLOCATE (lM%IEN(lM%eNoN,lM%nEl))
      DO e=1, lM%nEl
         READ(fid,*,END=001) i, lM%IEN(:,e)
      END DO
      lM%IEN = lM%IEN + 1
      CALL f%close()

      io%d = " Read a mesh with "//lM%nNo//" nodes and "//lM%nEl//
     2   " elements"

!     Satisfying right-handed property for the internal elements and
!     making sure everything is compatible
      CALL lM%check(.FALSE.)

      lM%ub = MAXVAL(lM%x,2)
      lM%lb = MINVAL(lM%x,2)

      ALLOCATE(rtmp(lM%nNo))
      rtmp = 1D0
      lM%vol = lM%integ(rtmp)
      DEALLOCATE(rtmp)

      io%o = " Total number of nodes: "//lM%nNo
      io%o = " Total number of elements: "//lM%nEl
      io%o = " Mesh lower bound: "//CLR(STR(lM%lb),3)
      io%o = " Mesh upper bound: "//CLR(STR(lM%ub),3)
      IF (nsd.EQ.2) THEN
         io%o = " Domain area: "//lM%vol
      ELSE
         io%o = " Domain volume: "//lM%vol
      END IF

      RETURN
 001  io%e = "A block of data is missing in the vtk file"
      END SUBROUTINE readPMsh
!---------------------------------------------------------------------
!     Check the mesh. If this is flag=nonL=.true., then we only 
!     check first 4 nodes of IEN array
      SUBROUTINE checkPMsh(lM, flag)
      CLASS(pMshType), INTENT(INOUT) :: lM
      LOGICAL, INTENT(IN) :: flag

      INTEGER Ac, b, i, sn(4), e, a, teNoN
      INTEGER, ALLOCATABLE :: incNodes(:)
      REAL(KIND=8), ALLOCATABLE :: v(:,:), xl(:,:)
       
      ALLOCATE(v(nsd,lM%eNoN), xl(nsd,lM%eNoN))
      
      teNoN = lM%eNoN
      IF (flag .AND. lM%eType.EQ.eType_BIQ) teNoN = 4
      ALLOCATE (incNodes(lM%nNo))
      incNodes = 0
      DO e=1, lM%nEl
         DO a=1, teNoN
            Ac = lM%IEN(a,e)
            IF (Ac.GT.lM%nNo .OR. Ac.LE.0) THEN
               io%e = "Element "//e//
     2            " contains an out of bound node"
            END IF
            incNodes(Ac) = 1
         END DO
      END DO
      DO Ac=1, lM%nNo
         IF (incNodes(Ac) .EQ. 0) THEN
            io%e = "Node "//Ac//
     2         " is isolated from the mesh"
         END IF
      END DO
      IF (lM%eType .EQ. eType_BRK) THEN
         io%o = "Make sure nodes in elements are arranged 1-4 on"
         io%o = "      one face and 5-8 on the opposite face"
         RETURN
      END IF

      DO e=1, lM%nEl
!     By default no change
         a = 1; b = 1
         IF (lM%eType.EQ.eType_BIL .OR. lM%eType.EQ.eType_BIQ) THEN
            xl     = lM%x(:,lM%IEN(1:4,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            v(:,3) = xl(:,4) - xl(:,3)
            v(:,4) = xl(:,1) - xl(:,4)
            sn(1)  = SGN(v(1,1)*v(2,2) - v(2,1)*v(1,2))
            sn(2)  = SGN(v(1,2)*v(2,3) - v(2,2)*v(1,3))
            sn(3)  = SGN(v(1,3)*v(2,4) - v(2,3)*v(1,4))
            sn(4)  = SGN(v(1,4)*v(2,1) - v(2,4)*v(1,1))
            i = sn(1) + sn(2) + sn(3) + sn(4)
            IF (i .EQ. 0) THEN
               IF (sn(1) .EQ. sn(2)) THEN
!     Sign is changing between 2-3 and 1-4
                  IF (sn(1) .EQ. 1) THEN
                     a = 1; b = 4
                  ELSE
                     a = 2; b = 3
                  END IF
               ELSE
!     Sign is changing between 1-2 and 3-4
                  IF (sn(1) .EQ. 1) THEN
                     a = 3; b = 4
                  ELSE
                     a = 1; b = 2
                  END IF
               END IF
            ELSE IF (i .EQ. -4) THEN
!     Sign is not changing, but this is going C.W.
               a = 1; b = 3
            ELSE IF (i.EQ.2 .OR. ANY(sn.EQ.0)) THEN
!     Two or more edges are on a straight line
               io%e = "Element "//e//" is distorted"
            END IF
         ELSE IF (lM%eType .EQ. eType_WDG) THEN
            xl     = lM%x(:,lM%IEN(:,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            v(:,3) = xl(:,4) - xl(:,1)
            v(:,4) = CROSS(v(:,1:2))
            sn(1)  = SGN(SUM(v(:,3)*v(:,4)))
            i = sn(1)
            IF (i .EQ. -1) THEN
!     Two nodes must be switched
               a = 1; b = 2
               Ac = lM%IEN(a,e)
               lM%IEN(a,e) = lM%IEN(b,e)
               lM%IEN(b,e) = Ac
               a = 4; b = 5
            ELSE IF (i .EQ. 0) THEN
!     Two or more edges are on a straight line
               io%e = "Element "//e//" is distorted"
            END IF
         ELSE IF (lM%eType .EQ. eType_TET) THEN
            xl     = lM%x(:,lM%IEN(:,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            v(:,3) = xl(:,4) - xl(:,3)
            v(:,4) = CROSS(v(:,1:2))
            sn(1)  = SGN(SUM(v(:,3)*v(:,4)))
            i = sn(1)
            IF (i .EQ. 1) THEN
!     Two nodes must be switched
               a = 1; b = 2
            ELSE IF (i .EQ. 0) THEN
!     Two or more edges are on a straight line
               io%e = "Element "//e//" is distorted"
            END IF
         ELSE IF (lM%eType .EQ. eType_TRI) THEN
            xl(:,1:3) = lM%x(:,lM%IEN(:,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            sn(1)  = SGN(v(1,1)*v(2,2) - v(2,1)*v(1,2))
            i = sn(1)
            IF (i .EQ. -1) THEN
!     Two nodes must be switched
               a = 1; b = 2
            ELSE IF (i .EQ. 0) THEN
!     Two or more edges are on a straight line
               io%e = "Element "//e//" is distorted"
            END IF
         END IF

         Ac = lM%IEN(a,e)
         lM%IEN(a,e) = lM%IEN(b,e)
         lM%IEN(b,e) = Ac
      END DO

      RETURN
      END SUBROUTINE checkPMsh
!---------------------------------------------------------------------
!     This routine integrate an equation over a particular domain
      FUNCTION integGPMsh(lM, s, lb, ub) RESULT(integ)
      CLASS(pMshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: s(:,:)
      INTEGER, INTENT(IN), OPTIONAL :: lb, ub
      REAL(KIND=8) integ

      INTEGER a, e, g, Ac, eNoN, fub, flb
      REAL(KIND=8) Jac, sHat, rtmp
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), Nx(:,:)
      
      IF (PRESENT(lb) .AND. PRESENT(ub)) THEN
         flb = lb
         fub = ub
      ELSE IF (.NOT.PRESENT(lb) .AND. PRESENT(ub)) THEN
         flb = ub
         fub = ub
      ELSE IF (PRESENT(lb) .AND. .NOT.PRESENT(ub)) THEN
         flb = lb
         fub = lb
      ELSE
         flb = 1
         fub = SIZE(s(:,1))
      END IF

      integ = 0D0
      eNoN  = lM%eNoN
      ALLOCATE(xl(nsd,eNoN), Nx(nsd,eNoN))
      DO e=1, lM%nEl
!     Updating the shape functions, if this is a NURB      
c         IF (msh(iM)%eType .EQ. eType_NRB) CALL NRBNNX(msh(iM), e)
         DO a=1, eNoN
            Ac      = lM%IEN(a,e)
            xl(:,a) = lM%x(:,Ac)
c            IF (mvMsh) xl(:,a) = xl(:,a) + Do(nsd+2:2*nsd+1,Ac)
         END DO

         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) CALL lM%dNdx(g,xl,Nx,Jac)
            sHat = 0D0
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               IF (flb .EQ. fub) THEN
                  rtmp = s(flb,Ac)
               ELSE
                  rtmp = NORM2(s(flb:fub,Ac))
               END IF
               sHat = sHat + rtmp*lM%N(a,g)
            END DO
            integ = integ + lM%w(g)*Jac*sHat
         END DO
      END DO

      RETURN
      END FUNCTION integGPMsh
!---------------------------------------------------------------------
!     This routine integrate an equation over a particular domain
      FUNCTION integSPMsh(lM, s) RESULT(integ)
      CLASS(pMshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN), TARGET :: s(:)
      REAL(KIND=8) integ

      REAL(KIND=8), POINTER :: sV(:,:)

      sV(1:1,1:SIZE(s)) => s

      integ = integGPMsh(lM, sV, 1)

      RETURN
      END FUNCTION integSPMsh
!---------------------------------------------------------------------
      SUBROUTINE freePMsh(this)
      CLASS(pMshType) :: this

      INTEGER i

      IF (ALLOCATED(this%fa)) THEN
         DO i=1, this%nFa
            CALL this%fa(i)%free()
         END DO
         DEALLOCATE(this%fa)
      END IF

      this%eType = eType_NA
      this%nEl   = 0
      this%nFa   = 0
      this%nNo   = 0

      IF (ALLOCATED(this%IEN)) DEALLOCATE(this%IEN)
      IF (ALLOCATED(this%lb))  DEALLOCATE(this%lb)
      IF (ALLOCATED(this%ub))  DEALLOCATE(this%ub)
      IF (ALLOCATED(this%x))   DEALLOCATE(this%x)

      RETURN 
      END SUBROUTINE freePMsh
      END MODULE PMSHMOD
