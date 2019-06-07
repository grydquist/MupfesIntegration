      MODULE PRTMOD
      USE INSMOD
      IMPLICIT NONE

!     Is used to initialize equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 2,
     2   itgO = 1,
     3   lsAl = 3,!LS_TYPE_PRT,
     5   sym  = "PRT")

!     Temporary holder for types of faces
      INTEGER, ALLOCATABLE :: faTyp(:)

!     Faces and their elements
      TYPE facelsType
      ! Elements of face
            INTEGER, ALLOCATABLE :: els(:)
      END TYPE

!     Individual box (elements)
      TYPE boxelType
!     Box dimensions (minx,maxx,miny,maxy,minz,maxz)
            REAL(KIND=8) :: dim(6)
!     Elements contained in searchbox
            INTEGER, ALLOCATABLE :: els(:)
!     Face elements in searchbox (face and element)
            TYPE(facelsType), ALLOCATABLE :: fa(:)
      END TYPE

!     Individual box (particles)
      TYPE boxpType
!     Box dimensions (minx,maxx,miny,maxy,minz,maxz)
            REAL(KIND=8) :: dim(6)
!     Total particles in box
            INTEGER :: nprt = 0
!     IDs of particles in box
            INTEGER, ALLOCATABLE :: c(:)
      END TYPE

!     Search boxes for elements
      TYPE sbeType
!     Created?
         LOGICAL :: crtd = .FALSE.
!     Number of boxes in each direction
         INTEGER n(3)
!     Size of boxes
         REAL(KIND=8) :: step(3)
!     Individual boxes
         TYPE(boxelType), ALLOCATABLE :: box(:)
!     Max and min sb values
         REAL(KIND=8) :: minx(3), maxx(3)
      CONTAINS
!     Sets up the search boxes pertaining to msh
         PROCEDURE :: new => newSbe
!     Returns serach box ID, provided the position of a point
         PROCEDURE :: id => idSbe
      END TYPE

!     Search boxes for particles
      TYPE sbpType
!     Created?
         LOGICAL :: crtd = .FALSE.
!     Number of boxes in each direction
         INTEGER n(3)
!     Total boxes
         INTEGER nt
!     Size of boxes
         REAL(KIND=8) :: step(3)
!     Individual boxes
         TYPE(boxpType), ALLOCATABLE :: box(:)
!     Max and min sb values
         REAL(KIND=8) :: minx(3), maxx(3)
      CONTAINS
!     Sets up the search boxes pertaining to msh
         PROCEDURE :: new => newSbp
!     Returns serach box ID, provided the position of a point
         PROCEDURE :: id => idSbp
!     Deallocates SB
         PROCEDURE :: free => freeSBp
!     Adds given particle IDs to SBs
         PROCEDURE :: addprt => addprtSBp
!     Removes given particle IDs to SBs
         PROCEDURE :: rmprt => rmprtSBp
      END TYPE

!     Your data type that you need to store in a Dynamic Sized Container
      TYPE pRawType
!     Eulerian element ID that this particle belongs to
         INTEGER(KIND=8) :: eID=1
!     Previous element ID
         INTEGER(KIND=8) :: eIDo=1
!     Particle ID
         INTEGER(KIND=8) :: pID=0
!     Searchbox ID (particles)
         INTEGER(KIND=8) :: sbIDp(8) = 0
!     Previous sbID (particles)
         INTEGER(KIND=8) :: sbIDpo(8) = 0
!     Searchbox ID (elements)
         INTEGER(KIND=8) :: sbIDe = 0
!     Previous sbID (elements)
         INTEGER(KIND=8) :: sbIDeo = 0
!     Crossed face ID
         INTEGER(KIND=8) :: faID(2) = 0
!     Position
         REAL(KIND=8) x(3)
!     Velocity
         REAL(KIND=8) u(3)
!     Previous Position
         REAL(KIND=8) :: xo(3) = 0D0
!     Previous Velocity
         REAL(KIND=8) :: uo(3) = 0D0
!     Position at collision
         REAL(KIND=8) :: xc(3) = 0D0
!     Velocity at collision
         REAL(KIND=8) :: uc(3) = 0D0
!     Shape functions in current element
         REAL(KIND=8), ALLOCATABLE :: N(:)
!     Remaining time in timestep after collision
         REAL(KIND=8) :: remdt
!     Time until collision
         REAL(KIND=8) :: ti
!     Other particles this particle has collided with
         INTEGER, ALLOCATABLE :: OthColl(:)
      END TYPE pRawType

!     Collection of particles
      TYPE, EXTENDS(eqType) :: prtType
!     Maximum capacity of container
         INTEGER :: maxN = 1024
!     Current size
         INTEGER :: n = 0
!     Current timestep
         REAL(KIND=8) :: dt = 0D0
!     Search boxes to find element hosting particles
         TYPE(sbeType) :: SBe
!     Search boxes for collisions
         TYPE(sbpType) :: SBp
!     Array pointing to full and empty entries
         INTEGER, ALLOCATABLE :: ptr(:)
!     Individual particle data
         TYPE(pRawType), ALLOCATABLE :: dat(:)
!     Velocity, from NS
         TYPE(varType), POINTER :: Uns => NULL()
!     Pressure, from NS
         TYPE(varType), POINTER :: Pns => NULL()
!     Material properties, from NS
         TYPE(matType), POINTER :: mns => NULL()
!     Two way coupling force passed to fluid
         TYPE(varType), POINTER :: twc => NULL()
!     Weighted volume
         REAL(KIND=8), ALLOCATABLE :: wV(:)
!     Collision counter
         INTEGER :: collcnt
!     List of collision pairs
         INTEGER, ALLOCATABLE :: collpair(:,:)
!     Time of collision pairs
         REAL(KIND=8), ALLOCATABLE :: collt(:)

      CONTAINS
!     Sets up all structure
         PROCEDURE :: setup => setupPrt
!     Returns shape function values at the location of a particle
         PROCEDURE :: shapeF => shapeFPrt
!     Seed the domai with particles
         PROCEDURE :: seed => seedPrt
!     Finds accelereation on a particle from drag
         PROCEDURE :: drag => dragPrt
!     Advance all particles one time step
         PROCEDURE :: solve => solvePrt
!     Detects collisions
         PROCEDURE :: findcoll => findcollPrt
!     Enacts collisions
         PROCEDURE :: collide => collidePrt
!     Finds which wall element the partice collides with
         PROCEDURE :: findwl => findwlPrt
!     Enacts collisions with walls
         PROCEDURE :: wall => wallPrt
!     Advances particle given time step
         PROCEDURE :: adv => advPrt

!     Evaulate the governing equation
         PROCEDURE :: eval3 => eval3Prt
         PROCEDURE :: eval2 => eval2Prt
!     Evaulate the governing equation boundary contribution
         PROCEDURE :: bEval => bEvalPrt
      END TYPE prtType

      INTERFACE prtType
            PROCEDURE :: newPrt
      END INTERFACE prtType

      CONTAINS

!####################################################################
!---------------------------------------------------------------------
      SUBROUTINE eval3Prt(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(prtType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)


      END SUBROUTINE eval3Prt

!---------------------------------------------------------------------
      PURE SUBROUTINE eval2Prt(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(prtType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)

      ! this is just a placeholder for now

      END SUBROUTINE eval2Prt
!---------------------------------------------------------------------
      PURE SUBROUTINE bEvalPrt(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(prtType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)


      END SUBROUTINE bEvalPrt

!---------------------------------------------------------------------
      SUBROUTINE writePrt(s, fName)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN) :: s
      CHARACTER(LEN=*), INTENT(IN) :: fName
      CHARACTER(LEN=3), PARAMETER :: names(2) = (/'Tmp','dID'/)
      INTEGER i, ip, l, fid, tN, pN, n
      INTEGER(KIND=8) pos(2)
      CHARACTER(LEN=stdL) hdr
      TYPE(pRawType), POINTER :: p
      TYPE(prtType) prt, tPrt
      INTEGER, ALLOCATABLE :: lN(:)
      REAL(KIND=8), ALLOCATABLE :: tmpX(:,:), tmpV(:,:), tmpT(:), 
     2   tmpI(:)


!     Work in progress... try looking at the dns one
      RETURN
      END SUBROUTINE writePrt
!---------------------------------------------------------------------
      FUNCTION newPrt(dmn, lst, Uns,mns,twc,Pns) RESULT(eq)
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst
      CLASS(varType), INTENT(IN), TARGET :: Uns
      CLASS(varType), INTENT(IN), TARGET :: Pns
      CLASS(matType), INTENT(IN), TARGET :: mns
      CLASS(varType), INTENT(IN), OPTIONAL, TARGET :: twc
      TYPE(prtType) :: eq
      TYPE(lstType), POINTER :: lPt1,lPt2,lPBC
      INTEGER nFa,iFa, typ2
      CHARACTER(LEN=stdL) stmp,typ1

      CALL eq%new(eqSp, dmn, lst)
      nFa  = lst%srch("Add BC")
      ALLOCATE(faTyp(nFa))
      eq%Uns => Uns
      eq%mns => mns
      eq%Pns => Pns
      IF(PRESENT(twc)) eq%twc => twc
      DO iFa=1, nFa
            lPBC => lst%get(stmp,"Add BC",iFa)
            lPt1 => lPBC%get(typ1,"Face type")
            SELECT CASE (LOWER(TRIM(typ1)))
            CASE ("inlet")
               faTyp(iFa) = 1
               lPt2 =>lPBC%get(typ2,"Number")
               eq%n = typ2
            CASE("outlet")
               faTyp(iFa) = 2
            CASE("wall")
               faTyp(iFa) = 3
            CASE DEFAULT
            io%e = "Select inlet, outlet, or wall Face type"         
         END SELECT
      ENDDO


      RETURN
      END FUNCTION newPrt

!--------------------------------------------------------------------
      SUBROUTINE setupPrt(eq, var)
      CLASS(prtType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
      INTEGER i,a, Ac
      REAL(KIND=8), ALLOCATABLE :: volt(:)

      open(88,file='pos.txt') 

      ALLOCATE(volt(eq%dmn%msh(1)%eNon))
      ALLOCATE(eq%dat(eq%n), eq%ptr(eq%n),eq%wV(eq%dmn%msh(1)%nNo))

      eq%wV = 0D0
      DO i=1, eq%n
         eq%ptr(i) = i
         ALLOCATE(eq%dat(i)%N(eq%dmn%msh(1)%eNoN))
      ENDDO
      eq%mat  => FIND_MAT('Particle')
      eq%var(1) = gVarType(nsd,'PVelocity',eq%dmn)
      eq%var(2) = gVarType(1,'PPosition',eq%dmn)

!     Assuming everything is in first mesh for searchboxes
      CALL eq%sbe%new(eq%dmn%msh(1))
!     Getting volumes of influence for each node
      DO i = 1,eq%dmn%msh(1)%nEl
            volt = effvolel(eq%dmn%msh(1),i)
            DO a = 1,eq%dmn%msh(1)%eNoN
                  Ac = eq%dmn%msh(1)%IEN(a,i)
                  eq%wV(Ac) = eq%wV(Ac) + volt(a)
            ENDDO
      ENDDO

      CALL eq%seed()

      END SUBROUTINE setupPrt
!-------------------------------------------------------------------- Elements
      SUBROUTINE newSbe(sb,msh)
      CLASS(mshType), INTENT(IN) :: msh
      CLASS(sbeType), INTENT(INOUT):: sb
      INTEGER :: ii,jj,cnt2,kk, iSb, xsteps(nsd)
     2  , xstepst(nsd,2), iSBmin, SBt,cx,cy,cz
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd),elvert(nsd,2),xzerot(nsd,2)
     2 , order(nsd,2), s,xzero(nsd)
      LOGICAL :: orderl(nsd)

      IF (sb%crtd) RETURN

!     Domain ranges
      diff(1)=MAXVAL(msh%x(1,:))-MINVAL(msh%x(1,:))
      diff(2)=MAXVAL(msh%x(2,:))-MINVAL(msh%x(2,:))
      diff(3)=MAXVAL(msh%x(3,:))-MINVAL(msh%x(3,:))

!     Ordering diff by index
      order(1,1) = MINLOC(diff,1)
      orderl = .true.
      orderl(order(1,1)) = .false.
      order(2,1) = MINLOC(diff,1, orderl)
      orderl(order(2,1)) = .false.
      order(3,1) = MINLOC(diff,1, orderl)
      order(:,2) = diff(order(:,1))
      order(:,2) = order(:,2)/order(1,2)

!     Scaling to get approximately cubic SBs equal to much greater than # elements
      s = (10*msh%nEl/(order(2,2)*order(3,2)))**(1D0/3D0)
      
!     First n estimate
      DO ii = 1,nsd
            sb%n(order(ii,1)) = INT(s*order(ii,2))
      ENDDO
      
!     Size of sb
      sb%step = diff/sb%n

!     Dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%box(sb%n(1)*sb%n(2)*sb%n(3)))

!     These sequences are just for allocating sbdim
      ALLOCATE(seq1(sb%n(3)*sb%n(2)),seq2(sb%n(3)*sb%n(1))
     2   ,seq3(sb%n(2)*sb%n(1)))

      seq1=(/(ii, ii=0, sb%n(2)*sb%n(3)-1, 1)/)*sb%n(1)+1
      cnt2=0
      DO ii=1,sb%n(1)*sb%n(3)
            seq2(ii)=ii+cnt2*(sb%n(2)-1)*sb%n(1)
            if (MOD(ii,sb%n(1)).eq.0) cnt2=cnt2+1
      ENDDO
      seq3=(/(ii, ii=0, sb%n(1)*sb%n(2)-1, 1)/)+1

      ! Direction 1
      DO ii=1,sb%n(1)
         sb%box(seq1+ii-1)%dim(1) = MINVAL(msh%x(1,:))
     2       + sb%step(1)*(ii-1)
      ENDDO

      ! Direction 2
      DO ii=1,sb%n(2)
         sb%box(seq2+(ii-1)*sb%n(1))%dim(3) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)
      ENDDO

      ! Direction 3
      DO ii=1,sb%n(3)
         sb%box(seq3+(ii-1)*sb%n(1)*sb%n(2))%dim(5)=
     2   MINVAL(msh%x(3,:)) + sb%step(3)*(ii-1)
      ENDDO

      sb%box%dim(2) = sb%box%dim(1) + sb%step(1)
      sb%box%dim(4) = sb%box%dim(3) + sb%step(2)
      sb%box%dim(6) = sb%box%dim(5) + sb%step(3)
      sb%minx(1) = minval(sb%box(:)%dim(1))
      sb%minx(2) = minval(sb%box(:)%dim(3))
      sb%minx(3) = minval(sb%box(:)%dim(5))
      

      ! Finding the sbs the box vertices are in
      DO ii=1,msh%Nel
            elvert(1,1) = MINVAL(msh%x(1,msh%IEN(:,ii)))
            elvert(1,2) = MAXVAL(msh%x(1,msh%IEN(:,ii)))
            elvert(2,1) = MINVAL(msh%x(2,msh%IEN(:,ii)))
            elvert(2,2) = MAXVAL(msh%x(2,msh%IEN(:,ii)))

         IF (nsd.eq.3) THEN
            elvert(3,1) = MINVAL(msh%x(3,msh%IEN(:,ii)))
            elvert(3,2) = MAXVAL(msh%x(3,msh%IEN(:,ii)))
         END IF

!           Set domain back to zero
            xzerot(1,:) = elvert(1,:) - sb%minx(1)
            xzerot(2,:) = elvert(2,:) - sb%minx(2)
            xzerot(3,:) = elvert(3,:) - sb%minx(3)

!           Number of searchbox steps in x,y,and z for both extreme vertices
            xstepst(1,:) = FLOOR(xzerot(1,:)/sb%step(1))
            xstepst(2,:) = FLOOR(xzerot(2,:)/sb%step(2))
            xstepst(3,:) = FLOOR(xzerot(3,:)/sb%step(3))

!           Difference between SB steps for extreme vertices
            xsteps = xstepst(:,2) - xstepst(:,1)

!           Furthest back SB the element is in
            iSBmin = xstepst(1,1) + sb%n(1)*xstepst(2,1) +
     2     sb%n(1)*sb%n(2)*xstepst(3,1) + 1

!           Now with this range, we can find all the SBs the element is in
!           Total SBs this particle is in
            SBt = (xsteps(1)+1)*(xsteps(2)+1)*(xsteps(3)+1)
            cx = 0
            cy = 0
            cz = 0

!           Loop over all SBs this element is in and add them
            DO jj = 1,SBt
!                 First we need to find the current ID of the SB we're in
                  iSB = iSBmin  + cx + cy*sb%n(1)
     2                          + cz*sb%n(1)*sb%n(2)
                  cx = cx + 1
                  IF (cx .gt. xsteps(1)) THEN
                        cx = 0
                        cy = cy + 1
                  END IF
                  IF (cy .gt. xsteps(2)) THEN
                        cy = 0
                        cz = cz + 1
                  END IF

!                 Add element to sb (if sb exists)
                  IF ((iSb.gt.0)
     2    .and.(iSb.le.sb%n(1)*sb%n(2)*sb%n(3))) THEN

!                 If this isn't the first element going in the box
                  IF (ALLOCATED(sb%box(iSb)%els))THEN
                        sb%box(iSb)%els = [sb%box(iSb)%els,ii]
!                 If this is the first element in the box
                  ELSE
                        ALLOCATE(sb%box(iSb)%els(1))
                        sb%box(iSb)%els(1) = ii
                  END IF
                  END IF
            ENDDO
      ENDDO

!     Same process as above, but for faces
      DO ii = 1,msh%nFa
            DO jj = 1,msh%fa(ii)%nEl

                  elvert(1,1) = 
     2             MINVAL(msh%x(1,msh%fa(ii)%IEN(:,jj)))
                  elvert(1,2) = 
     2             MAXVAL(msh%x(1,msh%fa(ii)%IEN(:,jj)))
         
                  elvert(2,1) = 
     2             MINVAL(msh%x(2,msh%fa(ii)%IEN(:,jj)))
                  elvert(2,2) = 
     2             MAXVAL(msh%x(2,msh%fa(ii)%IEN(:,jj)))
         
                  IF (nsd.eq.3) THEN
                        elvert(3,1) = 
     2             MINVAL(msh%x(3,msh%fa(ii)%IEN(:,jj)))
                        elvert(3,2) = 
     2             MAXVAL(msh%x(3,msh%fa(ii)%IEN(:,jj)))
                  END IF
         
!                 Set domain back to zero
                  xzerot(1,:) = elvert(1,:) - sb%minx(1)
                  xzerot(2,:) = elvert(2,:) - sb%minx(2)
                  xzerot(3,:) = elvert(3,:) - sb%minx(3)
         
!                 Find which searchbox the particle is in
!                 Number of searchbox steps in x,y,and z for both extreme vertices
                  xstepst(1,:) = FLOOR(xzerot(1,:)/sb%step(1))
                  xstepst(2,:) = FLOOR(xzerot(2,:)/sb%step(2))
                  xstepst(3,:) = FLOOR(xzerot(3,:)/sb%step(3))
         
!                 Difference between SB steps for extreme vertices
                  xsteps = xstepst(:,2) - xstepst(:,1)
         
!                 Furthest back SB the element is in
                  iSBmin = xstepst(1,1) + sb%n(1)*xstepst(2,1) +
     2                        sb%n(1)*sb%n(2)*xstepst(3,1) + 1
         
!                 Now with this range, we can find all the SBs the element is in
!                 Total SBs this particle is in
                  SBt = (xsteps(1)+1)*(xsteps(2)+1)*(xsteps(3)+1)
                  cx = 0
                  cy = 0
                  cz = 0

            DO kk = 1,SBt
                  iSB = iSBmin  + cx + cy*sb%n(1)
     2                          + cz*sb%n(1)*sb%n(2)
                  cx = cx + 1
                  IF (cx .gt. xsteps(1)) THEN
                        cx = 0
                        cy = cy + 1
                  END IF
                  IF (cy .gt. xsteps(2)) THEN
                        cy = 0
                        cz = cz + 1
                  END IF

!                 Add element to sb (if sb exists/element hasn't been added)
                  IF ((iSb.gt.0)
     2    .and.(iSb.le.sb%n(1)*sb%n(2)*sb%n(3))) THEN
!                 First we need to allocate the face structure of the sb...
                        IF(.not.ALLOCATED(sb%box(iSb)%fa))
     2                  ALLOCATE(sb%box(iSb)%fa(msh%nFa))             

!                 If this isn't the first element going in the box
                  IF (ALLOCATED(sb%box(iSb)%fa(ii)%els))THEN
                        sb%box(iSb)%fa(ii)%els =
     2                   [sb%box(iSb)%fa(ii)%els,jj]
!                 If this is the first element in the box
                  ELSE
                        ALLOCATE(sb%box(iSb)%fa(ii)%els(1))
                        sb%box(iSb)%fa(ii)%els = jj
                  END IF
                  END IF 
            ENDDO       
            ENDDO
      ENDDO

      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSbe

!-------------------------------------------------------------------- Particles
      SUBROUTINE newSbp(sb,msh,np,dmin)
      CLASS(mshType), INTENT(IN) :: msh
      INTEGER, INTENT(IN) :: np
      REAL(KIND=8), INTENT(IN) :: dmin(nsd)
      CLASS(sbpType), INTENT(INOUT):: sb
      INTEGER :: ii,cnt2, iSb(2**nsd),xsteps(2**nsd)
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), xzero(nsd), step(nsd), s

      IF (sb%crtd) RETURN
!     Smallest possible SB (diam + possible travel distance)

!     Domain ranges
      diff(1)=MAXVAL(msh%x(1,:))-MINVAL(msh%x(1,:))
      diff(2)=MAXVAL(msh%x(2,:))-MINVAL(msh%x(2,:))
      diff(3)=MAXVAL(msh%x(3,:))-MINVAL(msh%x(3,:))

!     Starting with particles the minimum traveling distance
      step = dmin
!     Scaling to get NB = NP
      s = (diff(1)/step(1)*diff(2)/step(2)*diff(3)/step(3)/np)
     2 **(1D0/3D0)

!     Applying this scaling
      step = step*MAX(s,1D0)
      
!     Getting the number of boxes from this
      sb%n = MAX(INT(diff/step),(/1,1,1/))
      
!     Size of sb
      sb%step = diff/((sb%n + 1D0)/2D0)
      sb%nt = sb%n(1)*sb%n(2)*sb%n(3)

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%box(sb%nt))

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(sb%n(3)*sb%n(2)),seq2(sb%n(3)*sb%n(1))
     2   ,seq3(sb%n(2)*sb%n(1)))

      seq1=(/(ii, ii=0, sb%n(2)*sb%n(3)-1, 1)/)*sb%n(1)+1
      cnt2=0
      DO ii=1,sb%n(1)*sb%n(3)
            seq2(ii)=ii+cnt2*(sb%n(2)-1)*sb%n(1)
            if (MOD(ii,sb%n(1)).eq.0) cnt2=cnt2+1
      ENDDO
      seq3=(/(ii, ii=0, sb%n(1)*sb%n(2)-1, 1)/)+1

      ! Allocating sb, such that they overlap by 50%
      ! Direction 1
      DO ii=1,sb%n(1)
         sb%box(seq1+ii-1)%dim(1) = MINVAL(msh%x(1,:)) 
     2       + sb%step(1)*(ii-1)/2
      ENDDO

      ! Direction 2
      DO ii=1,sb%n(2)
         sb%box(seq2+(ii-1)*sb%n(1))%dim(3) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)/2
      ENDDO

      ! Direction 3
      DO ii=1,sb%n(3)
         sb%box(seq3+(ii-1)*sb%n(1)*sb%n(2))%dim(5)=
     2   MINVAL(msh%x(3,:)) + sb%step(3)*(ii-1)/2
      ENDDO

      sb%box%dim(2) = sb%box%dim(1) + sb%step(1)
      sb%box%dim(4) = sb%box%dim(3) + sb%step(2)
      sb%box%dim(6) = sb%box%dim(5) + sb%step(3)
      sb%minx(1) = minval(sb%box(:)%dim(1))
      sb%minx(2) = minval(sb%box(:)%dim(3))
      sb%minx(3) = minval(sb%box(:)%dim(5))

!     Allocate SB particles in box
      DO ii = 1,sb%nt
            ALLOCATE(sb%box(ii)%c(1))
      ENDDO
      sb%box(:)%nprt = 0

      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSbp
!-------------------------------------------------------------------- Particles
!     Gets rid of SB
      SUBROUTINE freeSBp(sb)
      CLASS(sbpType), INTENT(INOUT):: sb
      
      IF (.not.(sb%crtd)) RETURN

      DEALLOCATE(sb%box)

      sb%crtd = .FALSE.


      END SUBROUTINE freeSBp
!-------------------------------------------------------------------- Particles
!     Put particles in SBs
      SUBROUTINE addprtSBp(sb,IDSBp,ip)
      CLASS(sbpType), INTENT(INOUT):: sb
      INTEGER, INTENT(IN) :: IDSBp(2**nsd), ip
      INTEGER :: i, cnted(2**nsd)

      cnted = 0
      DO i = 1,2**nsd
            IF ((IDSBp(i).le.sb%nt).and. (IDSBp(i).gt.0)
     2       .and. .not.(ANY(IDSBp(i).eq.cnted))) THEN
!     Increase number in box
                  sb%box(IDSBp(i))%nprt = sb%box(IDSBp(i))%nprt + 1

!     Add to end (a little different if it's the first element)
                  IF (sb%box(IDSBp(i))%nprt.ne.1) THEN
                        sb%box(IDSBp(i))%c =
     2                  [sb%box(IDSBp(i))%c,ip] 
                  ELSE
                        sb%box(IDSBp(i))%c(1) = ip
                  END IF
            END IF
            cnted(i) = IDSBp(i)
      ENDDO

      END SUBROUTINE addprtSBp
!-------------------------------------------------------------------- Particles
!     Put particles in SBs
      SUBROUTINE rmprtSBp(sb,IDSBp,ip)
      CLASS(sbpType), INTENT(INOUT):: sb
      INTEGER, INTENT(IN) :: IDSBp(2**nsd), ip
      INTEGER :: i, cnted(2**nsd)
      INTEGER, ALLOCATABLE :: tmpstck(:)
      
      cnted = 0
      DO i = 1,2**nsd
            IF ((IDSBp(i).le.sb%nt) .and. (IDSBp(i).gt.0)
     2       .and. .not.(ANY(cnted.eq.IDSBp(i)))) THEN
!     Decrease number in box
                  sb%box(IDSBp(i))%nprt = sb%box(IDSBp(i))%nprt - 1
!     Make array of box without particle
                  tmpstck = PACK(sb%box(IDSBp(i))%c, 
     2                           sb%box(IDSBp(i))%c.ne.ip)  
                  DEALLOCATE(sb%box(IDSBp(i))%c)
                  sb%box(IDSBp(i))%c = tmpstck
            END IF
            cnted(i) = IDSBp(i)
      ENDDO

      END SUBROUTINE rmprtSBp
!-------------------------------------------------------------------- Elements
!     Returns the ID of the searchboxes that contains point x
      FUNCTION idSBe(sb,x) RESULT(iSb)
      IMPLICIT NONE
      CLASS(sbeType), INTENT(IN) :: sb
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER iSb
      REAL(KIND=8) :: xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Set domain back to zero
      xzero(1) = x(1) - sb%minx(1)
      xzero(2) = x(2) - sb%minx(2)
      xzero(3) = x(3) - sb%minx(3)

      ! Find which searchbox the particle is in
      ! Number of searchbox steps in x,y,and z
      xsteps = FLOOR(xzero/sb%step)
      ! Check OOB
      IF (ANY(xsteps.ge.sb%n) .or. ANY(xsteps.lt.0)) THEN
            iSb = -1
            RETURN
      ENDIF
      ! SB it's in
      iSb = xsteps(1) + sb%n(1)*xsteps(2) +
     2   sb%n(1)*sb%n(2)*xsteps(3) + 1

      RETURN
      END FUNCTION idSBe
!--------------------------------------------------------------------
!     Returns the ID of the searchboxes that contains point x
      FUNCTION idSBp(sb,x) RESULT(iSb)
      IMPLICIT NONE
      CLASS(sbpType), INTENT(IN) :: sb
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER iSb(2**nsd)
      REAL(KIND=8) :: xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Set domain back to zero
      xzero(1) = x(1) - sb%minx(1)
      xzero(2) = x(2) - sb%minx(2)
      xzero(3) = x(3) - sb%minx(3)

      ! Find which searchbox the particle is in
      ! Number of searchbox steps in x,y,and z
      xsteps = FLOOR(2*xzero/sb%step)
      ! furthest searchbox in front
      iSb(1) = xsteps(1) + sb%n(1)*xsteps(2) +
     2   sb%n(1)*sb%n(2)*xsteps(3) + 1
      ! previous sb in x
      iSb(2) = iSb(1) - 1
      ! previous sb's in y
      iSb(3) = iSb(1) - sb%n(1)
      iSb(4) = iSb(3) - 1
      ! Next sb's in z (if available)
      if (nsd.eq.3) then
         iSb(5) = iSb(1) - sb%n(1)*sb%n(2)
         iSb(6) = iSb(5) - 1
         iSb(7) = iSb(5) - sb%n(1)
         iSb(8) = iSb(7) - 1
      end if

      RETURN
      END FUNCTION idSBp
!--------------------------------------------------------------------
!     Finds the element ID the particle is in. Also returns shape functions
      FUNCTION shapeFPrt(prt, ip, msh) RESULT(N)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN),TARGET :: prt
      INTEGER, INTENT(IN) :: ip
      TYPE(mshType), INTENT(IN) :: msh
      REAL(KIND=8) :: N(msh%eNoN),xl(nsd,msh%eNoN),Nt(msh%eNoN)
      TYPE(pRawType), POINTER :: p
      TYPE(boxelType),  POINTER :: b
      INTEGER :: ii, ind

      N= -1D0
      p => prt%dat(ip)
      IF (p%sbIDe .eq. -1) RETURN
      b => prt%sbe%box(p%sbIDe)
      DO ii=1,size(b%els)+1

            IF (ii.eq.1) THEN
                  ind = p%eIDo
            ELSE
                  ind = b%els(ii-1)
            END IF

            xl = msh%x(:,msh%IEN(:,ind))
            N = msh%nAtx(p%x,xl)

            ! Checking if all shape functions are positive
            IF (ALL(N.ge.-5D-7)) then
                  p%eID=ind
                  p%N = N
                  RETURN
            END IF
      ENDDO
         
      ! Catch to see if Sb is the issue
      DO ii = 1,msh%nEl
            xl = msh%x(:,msh%IEN(:,ii))
            Nt = msh%nAtx(p%x,xl)

            ! Checking if all shape functions are positive
            IF (ALL(Nt.ge.-5D-7)) then
                  print *, ii, ip, p%sbIDe,Nt
                  io%e = "SBe issue"
            END IF
      ENDDO

      ! If it loops through everything and Doesn't yield a positive shape function,
      ! the particle is outside the domain.
      RETURN
      END FUNCTION shapeFPrt
!--------------------------------------------------------------------
      SUBROUTINE seedPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      INTEGER ip
      TYPE(pRawType) p
      REAL, ALLOCATABLE :: N(:)
      REAL(KIND=8) dmin(nsd)

      ALLOCATE(N(prt%dmn%msh(1)%eNoN))
      dmin = prt%mat%D()
      CALL prt%sbp%new(prt%dmn%msh(1),prt%n,dmin)

      CALL RSEED(cm%id())
      DO ip=1, prt%n
           p%u    = 0D0
           p%pID  = INT(ip,8)
           CALL RANDOM_NUMBER(p%x(:))
           IF (nsd .EQ. 2) p%x(3) = 0D0
           p%x(3) = p%x(3)*15D0
           p%x = (p%x - (/0.5D0,0.5D0,0D0/))*20D0
           !p%x(1) = 0D0
           !p%x(2) = 0D0
           !p%x(3) = 150D0/ip
           !p%x(3) = 0.3D0
           !if (ip.eq.2) p%x(3)=0.1D0
           !if (mod(ip,2).eq.0) p%x(3) = 150/ip
           !if (mod(ip,2).eq.1) p%x(3) = 150/(ip+1)+.2D0
           p%sbIDp = prt%sbp%id(p%x)
           p%sbIDe = prt%sbe%id(p%x)
!          Reset other collisions
           IF(.not.ALLOCATED(p%OthColl)) ALLOCATE(p%OthColl(0))
           prt%dat(ip) = p
           N = prt%shapef(ip,prt%dmn%msh(1))
      ENDDO

      RETURN
      END SUBROUTINE seedPrt

!--------------------------------------------------------------------
!     Gets the acceleration due to drag on a single particle
      FUNCTION dragPrt(prt, ip) RESULT(apd)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER,INTENT(IN) :: ip
      TYPE(VarType),POINTER :: u
      TYPE(matType), POINTER :: mat
      TYPE(mshType), POINTER :: msh
      REAL(KIND=8) :: apd(nsd), fvel(nsd), taup, rhoP, mu,dp
      INTEGER :: ii, jj
      TYPE(pRawType), POINTER :: p
      ! Derived from flow
      REAL(KIND=8) :: fSN, magud, Rep, relvel(nsd), rhoF

      p => prt%dat(ip)
      dp = prt%mat%D()
      rhoP = prt%mat%rho()
      msh => prt%dmn%msh(1)
      u => prt%Uns

!     Fluid parameters
      rhoF = prt%mns%rho()
      mu  = prt%mns%mu()

!     Particle relaxation time
      taup = rhoP*dp**2D0/mu/18D0

!     Interpolate velocity at particle point
      fvel=0D0

      IF(.not. ALLOCATED(p%N))
     2 ALLOCATE(p%N(msh%eNoN)) 
      DO ii=1,nsd
         DO jj=1,msh%eNoN
            fvel(ii) = fvel(ii) + u%v(ii,msh%IEN(jj,p%eID))*p%N(jj)
         ENDDO
      ENDDO

      ! Relative velocity
      relvel = fvel - p%u
      ! Relative velocity magnitude
      magud = SUM(relvel**2D0)**0.5D0
      ! Reynolds Number
      Rep = dp*magud*rhoF/mu
      ! Schiller-Neumann (finite Re) correction
      fSN = 1D0 !+ 0.15D0*Rep**0.687D0   !!
      ! Stokes corrected drag force
      apD = fSN/taup*relvel

      END FUNCTION dragPrt
!--------------------------------------------------------------------
!     Detects and collisions and returns pairs and times
      SUBROUTINE findcollPrt(prt,id1,id2,lM)
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
      CLASS(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: id1, id2
      TYPE(pRawType), POINTER :: p1,p2

!     Calculating distance coefficient
      REAL(KIND=8) :: a, b, c, d, e, f, qa, qb,qc,zeros(2),tcr,dp
      REAL(KIND=8) :: Np1(nsd+1), Np2(nsd+1)
      INTEGER, ALLOCATABLE :: tmpcoll(:,:)

      p1 => prt%dat(id1)
      p2 => prt%dat(id2)
      dp  = prt%mat%D()

      p1%OthColl = [p1%OthColl,id2]
      p2%OthColl = [p2%OthColl,id1]

!     First, check if particles will collide at current trajectory
      a = p1%x(1) - p2%x(1)
      b = p1%u(1) - p2%u(1)
      c = p1%x(2) - p2%x(2)
      d = p1%u(2) - p2%u(2)
      if(nsd.eq.3) then
         e = p1%x(3) - p2%x(3)
         f = p1%u(3) - p2%u(3)
      else
         e=0D0
         f=0D0
      end if
      
      qa = b**2D0 + d**2D0 + f**2D0
      qb = 2D0*(a*b + c*d +e*f)
      qc = a**2D0 + c**2D0 + e**2D0 - ((dp + dp)/2D0)**2D0

!     Imaginary zeros means particles won't collide
      if ((qb**2D0-4D0*qa*qc).lt.0) RETURN

      ! Zeros are when the particle either enters or leaves vicinity of other particle
      zeros(1) = (-qb + sqrt(qb**2D0-4D0*qa*qc))/(2D0*qa)
      zeros(2) = (-qb - sqrt(qb**2D0-4D0*qa*qc))/(2D0*qa)

      ! Negative zeros mean the particle would collide previously in time
      if (ANY(zeros.le.0D0)) RETURN

      tcr = minval(zeros)

      ! Exit function if collision won't occur during (remaining)timestep
      if ((tcr.gt.p1%remdt) .and. (tcr.gt.p2%remdt)) 
     2 RETURN

      ! particle locations at point of collision
      p1%xc = p1%u*tcr + p1%x
      p2%xc = p2%u*tcr + p2%x
      p1%x = p1%xc ! I change these just so I can use shapeF
      p2%x = p2%xc ! I change them back below

!     Check if the particle is outside the domain
      p1%sbIDe = prt%sbe%id(p1%x)
      p2%sbIDe = prt%sbe%id(p2%x)     
      Np1 = prt%shapeF(id1, lM)
      Np2 = prt%shapeF(id2, lM)

      p1%x = p1%xo
      p2%x = p2%xo

!     OOB, no collisiion
      IF (ANY(Np1.lt.-5D-7) .or. ANY(Np2.lt.-5D-7)) THEN
            p1%sbIDe = prt%sbe%id(p1%x)
            p2%sbIDe = prt%sbe%id(p2%x)
            Np1 = prt%shapeF(id1, lM)
            Np2 = prt%shapeF(id2, lM)
            RETURN
      END IF

      p1%sbIDe = prt%sbe%id(p1%x)
      p2%sbIDe = prt%sbe%id(p2%x)
      Np1 = prt%shapeF(id1, lM)
      Np2 = prt%shapeF(id2, lM)

      p1%ti = tcr
      p2%ti = tcr                               !! just gets overwritten with mult collisions

      prt%collcnt = prt%collcnt+1
!     If collpair is allocated, this isn't the first collision to add
      IF (ALLOCATED(prt%collpair)) THEN
            ALLOCATE(tmpcoll(prt%collcnt,2))
            tmpcoll(1:prt%collcnt-1,1) = prt%collpair(1:prt%collcnt-1,1)
            tmpcoll(1:prt%collcnt-1,2) = prt%collpair(1:prt%collcnt-1,2)
            tmpcoll(prt%collcnt,1) = id1
            tmpcoll(prt%collcnt,2) = id2
            DEALLOCATE(prt%collpair)
            ALLOCATE(prt%collpair(prt%collcnt,2))
            prt%collpair = tmpcoll
            DEALLOCATE(tmpcoll)
!     IF it's not allocated, this is the first collision detected
      ELSE
            ALLOCATE(prt%collpair(prt%collcnt,2))
            prt%collpair(1,1) = id1
            prt%collpair(1,2) = id2
      END IF

!     Same procedure for collt
      IF (ALLOCATED(prt%collt)) THEN
            prt%collt = [prt%collt,tcr]! + (dt) !something like this  !! This needs to account for particles that have collided already
!     IF it's not allocated, this is the first collision detected
      ELSE
            ALLOCATE(prt%collt(prt%collcnt))
            prt%collt = tcr
      END IF

      RETURN
      END SUBROUTINE findcollPrt

!--------------------------------------------------------------------
      SUBROUTINE collidePrt(prt,id1,id2)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: id1, id2
      TYPE(pRawType), POINTER :: p1,p2
      REAL(KIND=8), ALLOCATABLE :: N(:)

!     Coefficients to make calculating parallel/perp vel easier
      REAL(KIND=8) :: vpar1, vpar2, vperp1, vperp2, dp, rho, k,mp
      REAL(KIND=8) :: n1(nsd), n2(nsd), t1(nsd), t2(nsd)
      REAL(KIND=8) :: pa, pb, temp(nsd)
      
      p1 => prt%dat(id1)
      p2 => prt%dat(id2)
      dp  = prt%mat%D()
      rho = prt%mat%rho()
      k   = prt%mat%krest()
      mp = pi*rho/6D0*dp**3D0

      ! Update location to collision location
      p1%x = p1%xc
      p2%x = p2%xc

      ! Vector parallel and pependicular to collision tangent line
      n1 = (p1%xc - p2%xc)/((dp + dp)/2)
      n2 = -n1
      temp = cross2(n1,p1%u)
      t1 = cross2(temp,n1)
      temp = cross2(n2,p2%u)
      t2 = cross2(temp,n2)

      ! Rare case with no perpendicular velocity
      if (ANY(ISNAN(t1))) t1 = 0D0
      if (ANY(ISNAN(t2))) t2 = 0D0
      
      ! Get precollision parallel and perpendicular velocities
      vperp1 = sum(t1*p1%u)
      vpar1  = sum(n1*p1%u)
      vperp2 = sum(t2*p2%u)
      vpar2  = sum(n2*p2%u)

      ! Note that perpendicular velocities DOn't change, so we only need to calculate parallel
      pa = mp*vpar1 - mp*vpar2
      pb = (-vpar1 - vpar2)*k

      vpar2 = (pa - mp*pb)/(mp + mp)
      vpar1 = pb + vpar2
      vpar2 = -vpar2

      ! V here is split into just two velocities, so just add them as vector

      p1%uc = vpar1*n1 + vperp1*t1
      p2%uc = vpar2*n2 + vperp2*t2

      p1%u = p1%uc
      p2%u = p2%uc

      p1%remdt = p1%remdt-p1%ti
      p2%remdt = p2%remdt-p2%ti

!     Update SB/element/shpfn info info of particles
      ALLOCATE(N(prt%dmn%msh(1)%eNoN))
      p1%sbIDp = prt%sbp%id(p1%x)
      p2%sbIDp = prt%sbp%id(p2%x)
      p1%sbIDe = prt%sbe%id(p1%x)
      p2%sbIDe = prt%sbe%id(p2%x)
      N = prt%shapeF(id1, prt%dmn%msh(1))
      N = prt%shapeF(id2, prt%dmn%msh(1))

      END SUBROUTINE collidePrt
!--------------------------------------------------------------------
      SUBROUTINE findwlPrt(prt,idp,msh)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      CLASS(mshType), INTENT(IN) :: msh
      TYPE(pRawType), POINTER :: p
      TYPE(boxelType),  POINTER :: b
      REAL(KIND=8) :: Jac, xXi(nsd,nsd), Am(nsd,nsd), x1(nsd), tc
      REAL(KIND=8) :: N(msh%eNoN),xi(nsd),Bm(nsd), xc(nsd) 
      INTEGER :: ii, jj, a,gEl, faceNS(1,msh%fa(1)%eNoN), cnt,kk,ll
     2 , faceN(msh%fa(1)%eNoN),facev
      REAL(KIND=8) s, t, mx, my, ux, uy, uz, lx, ly, lz, iJac,
     2 xl(nsd,msh%eNoN), magud

      p => prt%dat(idp)
      b => prt%sbe%box(p%sbIDe)

      p%faID = 0
      magud = SUM(p%u**2D0)**0.5D0

      faceloop: DO ii=1,msh%nFa
      DO jj=1,size(b%fa(ii)%els)
!           First we need to find the volumetric element associated with the face, and get x's associated with that
            gEl = msh%fa(ii)%gE(b%fa(ii)%els(jj))
            xl = msh%x(:,msh%IEN(:,gEl))
!           Next, we find which face of the volumetric element is the face element
            out: DO  kk= 1, msh%fa(ii)%eNoN
                  cnt = 1
                  DO ll = 1,msh%eNoN 
                        IF (msh%IEN(ll,gEl) .eq.
     2     msh%fa(ii)%IEN(kk,b%fa(ii)%els(jj))) THEN
                              faceNS(1,kk) = cnt
                              CYCLE out
                        ENDIF
                        cnt = cnt+1
                  ENDDO
            ENDDO out

            CALL QSORT(faceNS)
            faceN = faceNS(1,:)

            xXi = 0D0
            Bm = 0D0
            x1 = 0D0

!           Same process as NAtxEle for vol element
            DO a=1, msh%eNoN 
                  xXi(:,1) = xXi(:,1) +
     2          xl(:,a) * msh%Nx(1,a,1)
                  xXi(:,2) = xXi(:,2) +
     2          xl(:,a) * msh%Nx(2,a,1)
                  xXi(:,3) = xXi(:,3) +
     2          xl(:,a) * msh%Nx(3,a,1)   
!           Location of Gauss point (for Bm)
                  x1 = x1 + msh%N(a,1)*xl(:,a)
            ENDDO

            Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2      + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3      + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4      - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5      - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6      - xXi(1,3)*xXi(2,2)*xXi(3,1)
            iJac = 1D0/Jac

          Am(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))*iJac
          Am(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))*iJac
          Am(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))*iJac
          Am(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))*iJac
          Am(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))*iJac
          Am(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))*iJac
          Am(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))*iJac
          Am(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))*iJac
          Am(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))*iJac

!           Finding A*x_1
             DO a = 1,nsd
            Bm(a) = Am(a,1)*x1(1) + Am(a,2)*x1(2) + Am(a,3)*x1(3)
             ENDDO

!           Finding B = xi_1 - A*x_1
            Bm = msh%xi(:,1) - Bm

!     Now, knowing the face, we have information about the value of one parent coordinate
!     Which we can use to solve for tc. Then we can find the particle at point of collision 
!     with the plane in which the wall face is in. This is the same concept as NAtx, except
!     we're imposing the condition that xc will be on the same plane as the face we're checking

!     3D elements
      SELECT CASE(msh%eType)
      CASE(eType_BRK)
!     +y
      IF (ALL(faceN .eq. (/1,2,3,4/))) facev=1
!     +z
      IF (ALL(faceN .eq. (/1,2,5,6/))) facev=2
!     +x
      IF (ALL(faceN .eq. (/1,4,5,8/))) facev=3
!     -y
      IF (ALL(faceN .eq. (/5,6,7,8/))) facev=4
!     -z
      IF (ALL(faceN .eq. (/3,4,7,8/))) facev=5
!     -x
      IF (ALL(faceN .eq. (/2,3,6,7/))) facev=6

      SELECT CASE(facev)

      !     +y
            CASE(1)
      xi(2) = 1
      tc = -(-xi(2)+Am(2,1)*p%x(1) + Am(2,2)*p%x(2) + Am(2,3)*p%x(3)
     2  + Bm(2))/(Am(2,1)*p%u(1) + Am(2,2)*p%u(2) + Am(2,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)

      !     +z
            CASE(2)
      xi(3) = 1
      tc = -(-xi(3)+Am(3,1)*p%x(1) + Am(3,2)*p%x(2) + Am(3,3)*p%x(3)
     2  + Bm(3))/(Am(3,1)*p%u(1) + Am(3,2)*p%u(2) + Am(3,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

      !     +x
            CASE(3)
      xi(1) = 1
      tc = -(-xi(1)+Am(1,1)*p%x(1) + Am(1,2)*p%x(2) + Am(1,3)*p%x(3)
     2  + Bm(1))/(Am(1,1)*p%u(1) + Am(1,2)*p%u(2) + Am(1,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

      !     -y
            CASE(4)
      xi(2) = -1
      tc = -(-xi(2)+Am(2,1)*p%x(1) + Am(2,2)*p%x(2) + Am(2,3)*p%x(3)
     2  + Bm(2))/(Am(2,1)*p%u(1) + Am(2,2)*p%u(2) + Am(2,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)

      !     -z
            CASE(5)
      xi(3) = -1
      tc = -(-xi(3)+Am(3,1)*p%x(1) + Am(3,2)*p%x(2) + Am(3,3)*p%x(3)
     2  + Bm(3))/(Am(3,1)*p%u(1) + Am(3,2)*p%u(2) + Am(3,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

      !     -x
            CASE(6)
      xi(1) = -1
      tc = -(-xi(1)+Am(1,1)*p%x(1) + Am(1,2)*p%x(2) + Am(1,3)*p%x(3)
     2  + Bm(1))/(Am(1,1)*p%u(1) + Am(1,2)*p%u(2) + Am(1,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

      END SELECT

            ux = 1D0 + xi(1)
            uy = 1D0 + xi(2)
            uz = 1D0 + xi(3)
            lx = 1D0 - xi(1)
            ly = 1D0 - xi(2)
            lz = 1D0 - xi(3)
      
            N(1) = ux*uy*uz/8D0
            N(2) = lx*uy*uz/8D0
            N(3) = lx*uy*lz/8D0
            N(4) = ux*uy*lz/8D0
            N(5) = ux*ly*uz/8D0
            N(6) = lx*ly*uz/8D0
            N(7) = lx*ly*lz/8D0
            N(8) = ux*ly*lz/8D0
      CASE(eType_TET)
            IF (ALL(faceN .eq. (/1,2,3/))) facev = 1
            IF (ALL(faceN .eq. (/1,2,4/))) facev = 2
            IF (ALL(faceN .eq. (/1,3,4/))) facev = 3
            IF (ALL(faceN .eq. (/2,3,4/))) facev = 4

      SELECT CASE(facev)
            CASE(1)
!           This is the hard one. All xi sum to 1
            tc = (1 - p%x(1) * (Am(1,1)+Am(2,1)+Am(3,1))-
     2                p%x(2) * (Am(1,2)+Am(2,2)+Am(3,2))-
     3                p%x(3) * (Am(1,3)+Am(2,3)+Am(3,3))- 
     4                (Bm(1)+Bm(2)+Bm(3)))/
     5              ( p%u(1) * (Am(1,1)+Am(2,1)+Am(3,1))+
     6                p%u(2) * (Am(1,2)+Am(2,2)+Am(3,2))+
     7                p%u(3) * (Am(1,3)+Am(2,3)+Am(3,3)))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)
          
            CASE(2)
            xi(3) = 0
      tc = -(-xi(3)+Am(3,1)*p%x(1) + Am(3,2)*p%x(2) + Am(3,3)*p%x(3)
     2  + Bm(3))/(Am(3,1)*p%u(1) + Am(3,2)*p%u(2) + Am(3,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)          

            CASE(3)
            xi(2) = 0
      tc = -(-xi(2)+Am(2,1)*p%x(1) + Am(2,2)*p%x(2) + Am(2,3)*p%x(3)
     2  + Bm(2))/(Am(2,1)*p%u(1) + Am(2,2)*p%u(2) + Am(2,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)

            CASE(4)
            xi(1) = 0
      tc = -(-xi(1)+Am(1,1)*p%x(1) + Am(1,2)*p%x(2) + Am(1,3)*p%x(3)
     2  + Bm(1))/(Am(1,1)*p%u(1) + Am(1,2)*p%u(2) + Am(1,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)  

      END SELECT

            N(1) = xi(1)
            N(2) = xi(2)
            N(3) = xi(3)
            N(4) = 1D0 - xi(1) - xi(2) - xi(3)
      CASE(eType_WDG)
      io%e ="partcles not added for this element type"
            ux = xi(1)
            uy = xi(2)
            uz = 1D0 - xi(1) - xi(2)
            s  = (1D0 + xi(3))/2D0
            t  = (1D0 - xi(3))/2D0
                  
            N(1) = ux*t
            N(2) = uy*t
            N(3) = uz*t
            N(4) = ux*s
            N(5) = uy*s
            N(6) = uz*s      
!     2D elements         
      CASE(eType_TRI)
      io%e="partcles not added for this element type"
         N(1) = xi(1)
         N(2) = xi(2)
               N(3) = 1D0 - xi(1) - xi(2)
      CASE(eType_BIL)
      io%e="partcles not added for this element type"
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
               
         N(1) = ux*uy/4D0
         N(2) = lx*uy/4D0
         N(3) = lx*ly/4D0
         N(4) = ux*ly/4D0
      CASE(eType_BIQ)
      io%e="partcles not added for this element type"
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         mx = xi(1)
         my = xi(2)
               
         N(1) =  mx*lx*my*ly/4D0
         N(2) = -mx*ux*my*ly/4D0
         N(3) =  mx*ux*my*uy/4D0
         N(4) = -mx*lx*my*uy/4D0
         N(5) = -lx*ux*my*ly/2D0
         N(6) =  mx*ux*ly*uy/2D0
         N(7) =  lx*ux*my*uy/2D0
         N(8) = -mx*lx*ly*uy/2D0
         N(9) =  lx*ux*ly*uy
      
!     1D elements         
      CASE(eType_LIN)
      io%e="partcles not added for this element type"
         N(1) = (1D0 - xi(1))/2D0
         N(2) = (1D0 + xi(1))/2D0
      CASE(eType_QUD)
      io%e="partcles not added for this element type"
         N(1) = -xi(1)*(1D0 - xi(1))/2D0
         N(2) =  xi(1)*(1D0 + xi(1))/2D0
         N(3) = (1D0 - xi(1))*(1D0 + xi(1))
      END SELECT

! End NatxiEle

            IF (ALL(N.ge.-5D-7).and. (tc.ge.-5D-7*magud)
     2      .and. tc.lt.prt%dt) THEN
                  p%faID(1) = ii
                  p%faID(2) = b%fa(ii)%els(jj)
                  p%ti = tc
                  RETURN
            ENDIF

      ENDDO
      ENDDO faceloop
      print *, p%N, p%x, p%u, p%sbIDe, idp
      io%e = "Wrong searchbox"
      
      END SUBROUTINE findwlPrt
!--------------------------------------------------------------------
      SUBROUTINE wallPrt(prt, idp, msh)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      CLASS(mshType), INTENT(IN) :: msh
      TYPE(pRawType), POINTER :: p
      INTEGER :: ii, rndi, jj, Ac
      REAL(KIND=8) :: dp, rho, k, nV(nsd), tV(nsd),
     2 a(nsd), b(nsd), vpar, vperp, temp(nsd), rnd,
     3 apd(nsd), mp, rhoF

      p   =>prt%dat(idp)
      dp  = prt%mat%D()
      rho = prt%mat%rho()
      k   = prt%mat%krest()
      rhoF  = prt%mns%rho()
      mp = pi*rho/6D0*dp**3D0

!     Get first order drag
      apd = prt%drag(idp)
      p%u = p%u + apd*p%ti

!     Send drag to fluid
      DO jj=1,msh%eNoN
            Ac = msh%IEN(jj,p%eID)
            prt%twc%v(:,Ac) = prt%twc%v(:,Ac) +
     2      apd*mP/rhoF/prt%wV(Ac)*p%N(jj)
      ENDDO

!     Advance to collision location
      p%xc = p%u*p%ti + p%x
      p%x  = p%xc
      p%remdt = p%remdt - p%ti

      IF (faTyp(p%faID(1)) .EQ. 3) THEN
!     Hitting wall

!     Get normal/tangent vector
            a = msh%x(:,msh%fa(p%faID(1))%IEN(1,p%faID(2))) - 
     2    msh%x(:,msh%fa(p%faID(1))%IEN(2,p%faID(2)))
            b = msh%x(:,msh%fa(p%faID(1))%IEN(1,p%faID(2))) - 
     2    msh%x(:,msh%fa(p%faID(1))%IEN(3,p%faID(2)))
            nV = CROSS2(a,b)
            temp = CROSS2(nV,p%u)
            tV = CROSS2(temp,nV)
      ! Rare case with no perpendicular velocity
            IF (ANY(ISNAN(tV))) tV = 0D0
            vperp = sum(tV*p%u)
            vpar  = sum(nV*p%u)*k
      !     Change velocities to after collision
            p%u = -vpar*nV +vperp*tV

      ELSE

!     Exiting domain
      faloop:DO ii = 1,msh%nFa
!           Search for inlet to put particle back into
            IF (faTyp(ii) .EQ. 1) THEN
!           Select random node on face to set as particle position
            CALL RANDOM_NUMBER(rnd)
            rndi = FLOOR(msh%fa(ii)%nEl*rnd + 1)
            p%x = msh%x(:,msh%fa(ii)%IEN(1,rndi))*.99       !! silly hack b/c particles on edges get lost, should be fixed w/ periodic

            EXIT faloop
            END IF
      ENDDO faloop

      END IF

      END SUBROUTINE wallPrt

!--------------------------------------------------------------------
      ! Finds the volume influenced by moving one node in one element,
      ! for all nodes in that element
      FUNCTION effvolel(msh,e) RESULT(effvol)
      IMPLICIT NONE
      TYPE(mshtype), INTENT(IN) :: msh
      INTEGER, INTENT(IN) :: e
      REAL(KIND=8) :: effvol(msh%eNoN)
      INTEGER g,a, Ac(msh%eNoN)
      REAL(KIND=8) :: Jac, x(nsd, msh%eNoN), Nx(nsd,msh%eNoN)

      effvol = 0D0

!     Indices and coordinates of nodes
      DO a = 1,msh%eNoN
            Ac(a) = msh%IEN(a,e)
      ENDDO

      x = msh%x(:,Ac)      

      DO g = 1, msh%nG
!     First, we want the Jacobian, which (if shpfns are linear) is the same for all gauss points
      IF (g.EQ.1 .OR. .NOT.msh%lShpF) CALL msh%dNdx(g,x,Nx,Jac)
            DO a = 1, msh%eNoN
!           Numerically integrate to get volume of each node
                  effvol(a) = effvol(a) + msh%N(a,g)*Jac*msh%w(g)
            ENDDO
      ENDDO

      END FUNCTION

!--------------------------------------------------------------------
      ! I use cross products a couple times above and didn't want to integrate the util one
      ! Also normalizes to unit vector
      FUNCTION cross2(v1,v2) RESULT(cross)
      IMPLICIT NONE
      REAL(KIND=8) :: cross(nsd)
      REAL(KIND=8), INTENT(IN) :: v1(nsd), v2(nsd)

      cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
      cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
      cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

      cross = cross/sqrt(cross(1)**2D0 + cross(2)**2D0 + cross(3)**2D0)
      
      END FUNCTION cross2
!--------------------------------------------------------------------
      SUBROUTINE advPrt(prt, idp)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      INTEGER, ALLOCATABLE :: tmpstck(:)
      TYPE(matType), POINTER :: mat
      TYPE(mshType), POINTER :: msh
      TYPE(pRawType), POINTER :: p, tp
      TYPE(prtType), TARGET :: tmpprt
      TYPE(sbpType), POINTER :: sbp
      TYPE(sbeType), POINTER :: sbe
      TYPE(VarType),POINTER :: u
      REAL(KIND=8) rhoF,
     2   mu, mp, g(nsd), dp,
     3   rhoP, Ntmp(nsd+1), tt(4)
      REAL(KIND=8), ALLOCATABLE :: N(:)
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd)
      REAL(KIND=8) :: prtxpred(nsd), pvelpred(nsd), tmpwr(nsd),mom
      INTEGER :: a, Ac, i, tsbIDp(2**nsd),tsbIDe, teID

!     Gravity
!! Find where grav actually is? (maybe mat%body forces)
      g=0D0
      g(3)=-1D0/idp
      !if (mod(idp,2).eq.0) g(3) =1D0
      !if (mod(idp,2).eq.1) g(3) =-1D0
      !if (idp .eq.2) g(3) = 1D0
      !tmpprt = prt                  !! Expensive with more sb
      msh => prt%dmn%msh(1)
      ALLOCATE(N(msh%eNoN))
      p  => prt%dat(idp)
      sbp => prt%sbp
      sbe => prt%sbe
      tp => tmpprt%dat(1)
      u => prt%Uns
      !tmpprt%sb = sb                !! Expensive with more sb
      rhoF  = prt%mns%rho()
      rhoP  = prt%mat%rho()
      mu    = prt%mns%mu()
      dp    = prt%mat%D()
      mp = pi*rhoP/6D0*dp**3D0

!     Get drag acceleration
      apd = prt%drag(idp)
!     Total acceleration (just drag and buoyancy now)
      apT = apd + g*(1D0 - rhoF/rhoP)

!!    For now, I'm just going to DO 1st order. It eases a lot of optimization stuff

!!     2nd order advance (Heun's Method)
!!     Predictor
!      pvelpred = p%u + p%remdt*apT
!      prtxpred = p%x + p%remdt*p%u
!      tp%u  = pvelpred
!      tp%x  = prtxpred
!
!!     Find which searchbox prediction is in
!      tp%sbID = sb%id(tp%x)   !! For some reason this changes apT???? Does not make any sense to me at all
!!     Get shape functions/element of prediction
!      Ntmp = tmpprt%shapeF(1, msh)
!!     Check if predictor OOB
!      IF (ANY(Ntmp.le.-1D-7)) THEN
!!     This should advance to the edge of the wall, and then change the velocitry, as well as giving remdt
!            CALL prt%findwl(idp,msh)
!            CALL prt%wall(idp,msh)
!            p%x = p%x + p%remdt*p%u
!            p%wall = .TRUE.
!            RETURN
!      END IF
!!     Total acceleration (just drag and buoyancy now) ! again, because above changes this for some reason !!
!      apT = apd + g*(1D0 - rhoF/rhoP)
!
!!     Get drag acceleration of predicted particle
!      apdpred = tmpprt%drag(1)
!      apTpred = apdpred + g*(1D0 - rhoF/rhoP)
!!     Corrector
      p%u = p%u + p%remdt*apT!0.5D0*p%remdt*(apT+apTpred)
      p%x = p%remdt*p%u + p%x

!     Send drag to fluid
      DO a=1,msh%eNoN
            Ac = msh%IEN(a,p%eID)
            prt%twc%v(:,Ac) = prt%twc%v(:,Ac) +
     2      apd*mP/rhoF/prt%wV(Ac)*p%N(a)
!     2      0.5D0*(apd + apdpred)*mP/rhoF/prt%wV(Ac)*p%N(a)
      ENDDO

      IF (prt%itr .EQ. 0) THEN
            tmpwr = apd
            write(88,*) sqrt(tmpwr(1)**2+tmpwr(2)**2+tmpwr(3)**2)*mp
            !print *, sqrt(tmpwr(1)**2+tmpwr(2)**2+tmpwr(3)**2)*mp
            !mom =prt%dmn%msh(1)%integ(u%v, 3)
            !print *, mom
            !call sleep(1)
      END IF

!     Check if particle went out of bounds
      tsbIDe = p%sbIDe
      tsbIDp = p%sbIDp
      teID = p%eID
      p%sbIDp = sbp%id(p%x)
      p%sbIDe = sbe%id(p%x)
      N = prt%shapeF(idp, msh)

!     If so, do a wall collision and continue on
      IF (ANY(N .le. -5D-7)) THEN
            p%x = p%xo
            p%sbIDe = tsbIDe
            p%sbIDp = tsbIDp
            p%eID = teID
            CALL prt%findwl(idp,msh)
            CALL prt%wall(idp,msh)
            p%x = p%x + p%remdt*p%u
            p%sbIDe = sbe%id(p%x)
            p%sbIDp = sbp%id(p%x)
            N = prt%shapeF(idp, msh)
      END IF

      RETURN

      END SUBROUTINE advPrt
!--------------------------------------------------------------------
      SUBROUTINE solvePrt(eq)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT):: eq
      INTEGER ip, a, e, eNoN, i ,j, k,l, subit, citer,i2,i1
     2 ,IDSBp(2**nsd),IDSBp1(2**nsd),IDSBp2(2**nsd), i12,i22
      TYPE(mshType), POINTER :: lM
      INTEGER, ALLOCATABLE :: tmpstck(:), N(:)
      REAL(KIND=8):: dtp,maxdtp,sbdt(nsd),dp,taup,rhoP,mu,tim,
     2 P1, P2, rhoF, umax(3) =0D0, dmin(nsd), tic, toc
      CHARACTER (LEN=stdl) fName

      lM => eq%dmn%msh(1)
      rhoP  = eq%mat%rho()
      rhoF  = eq%mns%rho()
      mu    = eq%mns%mu()
      dp    = eq%mat%D()
      ALLOCATE(N(lm%eNoN))

!     Reset twc force to zero
      eq%twc%v(:,:) = 0D0
!      eq%twc%v(3,1952) = 1D0

!     Reset collision counter
      eq%collcnt = 0

!     Particle relaxation time
      taup = rhoP*dp**2D0/mu/18D0

!     Appropriate time step
      dtp = MIN(taup,dt)
!     Subiterations
      subit = 1!FLOOR(dt/dtp)
      dtp = dt/subit
      eq%dt = dtp

      DO l = 1,subit

      IF (eq%itr .EQ. 0) THEN
!     Reset SB for particles
            DO i = 1, eq%n
                  umax = MAX(umax,ABS(eq%dat(i)%u))
            ENDDO
            dmin = umax*dtp+dp
            CALL eq%sbp%free()
            CALL eq%sbp%new(eq%dmn%msh(1), eq%n,dmin)
      END IF

      DO i=1,eq%n
!     If it's the first iteration, update last velocity, position, info
            IF (eq%itr .EQ. 0) THEN
                  eq%dat(i)%xo = eq%dat(i)%x
                  eq%dat(i)%uo = eq%dat(i)%u
                  eq%dat(i)%eIDo = eq%dat(i)%eID
                  eq%dat(i)%sbIDeo= eq%dat(i)%sbIDe
                  eq%dat(i)%sbIDpo= eq%dat(i)%sbIDp

!                 Put particles in SBs
                  idSBp = eq%dat(i)%sbIDp
                  CALL eq%sbp%addprt(idSBp,i)
            ENDIF  

!     Set position and velocity to old variables in preparation for iteration with INS
            eq%dat(i)%x = eq%dat(i)%xo
            eq%dat(i)%u = eq%dat(i)%uo
            eq%dat(i)%eID = eq%dat(i)%eIDo
            eq%dat(i)%sbIDe = eq%dat(i)%sbIDeo
            eq%dat(i)%sbIDp = eq%dat(i)%sbIDpo

!     Set initial advancing time step to solver time step
            eq%dat(i)%remdt = dtp
!     Reset collisions from last time step
            DEALLOCATE(eq%dat(i)%OthColl)
            ALLOCATE  (eq%dat(i)%OthColl(0))
      ENDDO

      tic = CPUT()
      DO i = 1,eq%n
            idSBp = eq%dat(i)%sbIDp
            DO j=1,2**nsd
            IF ((IDSBp(j).gt.0).and. (IDSBp(j).le.eq%sbp%nt)) THEN
                  DO k=1,eq%sbp%box(IDSBp(j))%nprt
                        i2 = eq%sbp%box(IDSBp(j))%c(k)        
!                       Check if the particle collides with any other particles during step. Then add to list if so
                        IF ((i.ne.i2) .and. 
     2                   .not.ANY(eq%dat(i)%OthColl.eq.i2)) THEN
                              CALL eq%findcoll(i,i2,lM)
                        ENDIF
                  ENDDO
            END IF
            ENDDO
      ENDDO
      toc = CPUT()
      toc = toc - tic
      print *, toc
      if (cm%mas()) then
         open(123,file='speed_'//STR(eq%n)//'.txt',position='append')
         write(123,*) toc
         close(123)
      end if 

      IF (ALLOCATED(eq%collt)) THEN
!     Keep looping over collisions, adding more as they appear, until there aren't any left
      DO WHILE(MINVAL(eq%collt).le.dtp)
!     Enact collisions, starting with shortest time to collision, check for more collisions, add any you find in
            citer = MINLOC(eq%collt, 1)
            i1 = eq%collpair(citer,1)
            i2 = eq%collpair(citer,2)

!           Remove particles from sb before collision
            idSBp = eq%dat(i1)%sbIDp
            CALL eq%sbp%rmprt(idSBp,i1)
            idSBp = eq%dat(i2)%sbIDp
            CALL eq%sbp%rmprt(idSBp,i2)

!           Enact collisions
            CALL eq%collide(i1,i2)

!           Put particles into their new SBs
            idSBp = eq%dat(i1)%sbIDp
            CALL eq%sbp%addprt(idSBp,i1)
            idSBp = eq%dat(i2)%sbIDp
            CALL eq%sbp%addprt(idSBp,i2)

            eq%collt(citer) = dtp + 1
!     Get rid of any previous collisions with these particles
            WHERE(eq%collpair(:,1).eq.i1) eq%collt = dtp + 1
            WHERE(eq%collpair(:,2).eq.i1) eq%collt = dtp + 1
            WHERE(eq%collpair(:,1).eq.i2) eq%collt = dtp + 1
            WHERE(eq%collpair(:,2).eq.i2) eq%collt = dtp + 1

            DEALLOCATE(eq%dat(i1)%OthColl)
            ALLOCATE(  eq%dat(i1)%OthColl(1))
            eq%dat(i1)%OthColl(1) = i2
            DEALLOCATE(eq%dat(i2)%OthColl)
            ALLOCATE(  eq%dat(i2)%OthColl(1))
            eq%dat(i2)%OthColl(1) = i1

            idSBp1 = eq%dat(i1)%sbIDp
            idSBp2 = eq%dat(i2)%sbIDp
            DO i = 1,2**nsd

            IF ((IDSBp1(i).gt.0).and. (IDSBp1(i).le.eq%sbp%nt)) THEN
                  DO k=1,eq%sbp%box(IDSBp1(i))%nprt
                        i12 = eq%sbp%box(IDSBp1(i))%c(k)        
!                       Check if the particle collides with any other particles after initial collision. Then add to list if so
                        IF ((i1.ne.i12) .and. 
     2                   .not. ANY(eq%dat(i1)%OthColl.eq.i12))
     3                   CALL eq%findcoll(i1,i12,lM)
                  ENDDO
            END IF

            IF ((IDSBp2(i).gt.0).and. (IDSBp2(i).le.eq%sbp%nt)) THEN
                  DO k=1,eq%sbp%box(IDSBp2(i))%nprt
                        i22 = eq%sbp%box(IDSBp2(i))%c(k)        
!                       Check if the particle collides with any other particles after initial collision. Then add to list if so
                        IF ((i2.ne.i22) .and.
     2                   .not. ANY(eq%dat(i2)%OthColl.eq.i22))
     3                   CALL eq%findcoll(i2,i22,lM)
                  ENDDO
            END IF
            ENDDO
      ENDDO
      ENDIF

!     Reset collpair/collt
      IF (ALLOCATED(eq%collpair)) DEALLOCATE(eq%collpair, eq%collt)

      tic = CPUT()
      DO i = 1,eq%n
!           Now no particles will collide on their path. Safe to advance them 
            CALL eq%adv(i)
!            IF (eq%itr .EQ. 0)  print *, eq%dat(i)%x(3),
!     2            eq%dat(i)%u(3)
      ENDDO
      toc = CPUT()
      toc = toc - tic
      print *, toc
            
      P1 = eq%dmn%msh(1)%integ(1,eq%Pns%s)
      P2 = eq%dmn%msh(1)%integ(2,eq%Pns%s)

      IF (eq%itr .EQ. 0) write(88,*) eq%dat(1)%u(3), !eq%dat(2)%u
     2 eq%dmn%msh(1)%integ(eq%Uns%v, 3)*rhoF,
     3 time
!     3  (eq%dmn%msh(1)%integ(1, eq%Uns%v,3 ) -
!     4  eq%dmn%msh(1)%integ(2, eq%Uns%v,3 ) -
!     5  eq%dmn%msh(1)%integ(3, eq%Uns%v,3 ))*1.2D0,
     6  ,P1,P2
      ENDDO

!     Write particle data if it's a multiple of ten timestep
      IF (mod(cTs,2).eq.0 .and. eq%itr.eq.0) THEN
            fName = "prt_"//STR(cTs)//".vtk"
            CALL writePrt(eq,fName)
      ENDIF

!     Updating norm for solution control
      CALL eq%upNorm(eq%ls%RI%iNorm)
      
!     Checking for exceptions
      CALL io%w%checkException()

      RETURN
      END SUBROUTINE solvePrt
      
      END MODULE PRTMOD


      !! Urgent fixes:
      !! Mult collisions still need quite a bit of work to keep time consistent (if remdt1 .ne. remdt2)
      !! Don't have it implemented so it can hit multiple walls
      !! Make sbs not order NB^1/3