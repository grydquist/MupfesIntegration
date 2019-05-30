      MODULE PRTMOD
      USE INSMOD
!     FSI is using PRTMOD
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
         LOGICAL, ALLOCATABLE :: OthColl(:)
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
!     Data
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
      END DO


      RETURN
      END FUNCTION newPrt

!--------------------------------------------------------------------
      SUBROUTINE setupPrt(eq, var)
      CLASS(prtType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
      INTEGER i,a, Ac,cnt
      REAL(KIND=8), ALLOCATABLE :: volt(:)

      open(88,file='pos.txt') 

      ALLOCATE(volt(eq%dmn%msh(1)%eNon))
      ALLOCATE(eq%dat(eq%n), eq%ptr(eq%n),eq%wV(eq%dmn%msh(1)%nNo))
      ALLOCATE(eq%collt(eq%n**2),eq%collpair(eq%n**2,2))

      eq%wV = 0D0
      cnt =0
      DO i=1, eq%n
         eq%ptr(i) = i
         ALLOCATE(eq%dat(i)%N(eq%dmn%msh(1)%eNoN))
      END DO
      eq%mat  => FIND_MAT('Particle')
      eq%var(1) = gVarType(nsd,'PVelocity',eq%dmn)
      eq%var(2) = gVarType(1,'PPosition',eq%dmn)

!     Assuming everything is in first mesh for searchboxes
      CALL eq%sbe%new(eq%dmn%msh(1), eq%n)
      CALL eq%sbp%new(eq%dmn%msh(1), eq%n)
!     Getting volumes of influence for each node
      DO i = 1,eq%dmn%msh(1)%nEl
            volt = effvolel(eq%dmn%msh(1),i)
            DO a = 1,eq%dmn%msh(1)%eNoN
                  Ac = eq%dmn%msh(1)%IEN(a,i)
                  eq%wV(Ac) = eq%wV(Ac) + volt(a)
            END DO
      END DO

      CALL eq%seed()

      END SUBROUTINE setupPrt
!-------------------------------------------------------------------- Elements
      SUBROUTINE newSbe(sb,msh,np)
      CLASS(mshType), INTENT(IN) :: msh
      INTEGER, INTENT(IN) :: np
      CLASS(sbeType), INTENT(INOUT):: sb
      TYPE(boxelType) :: tmpbox
      INTEGER :: ii,jj,cnt2,kk,ll, iSb, xsteps(nsd)
     2  , n(nsd)
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), elbox(2*nsd,msh%nEl),
     2 elboxf(2**nsd),elvert(nsd,2**nsd),xzero(nsd), order(nsd,2)
     3 , maxel(nsd), s
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
      order(2,1) = MINLOC(diff,1, MASK = orderl)
      orderl(order(2,1)) = .false.
      order(3,1) = MINLOC(diff,1, MASK = orderl)
      order(:,2) = diff(order(:,1))
      order(:,2) = order(:,2)/order(1,2)

!     Scaling to get approximately cubic SBs equal to approx number of particles
      s = ((np*msh%nEl)**(0.5)/(order(2,2)*order(3,2)))**(1D0/3D0)
      
!     First n estimate
      DO ii = 1,nsd
            n(order(ii,1)) = INT(s*order(ii,2))
      ENDDO
      
!     Size of sb
      sb%step = diff/n

!     Now we check to make sure the dimensions of the SBs are bigger than the elements, so we don't miss an element
!     We start by making boxes around the elements, which we will use later as well, and getting max element size
      
      maxel = 0
      DO ii=1,msh%Nel
            DO jj=1,nsd
               elbox(2*jj-1,ii) = MINVAL(msh%x(jj,msh%IEN(:,ii)))
               elbox(2*jj  ,ii) = MAXVAL(msh%x(jj,msh%IEN(:,ii)))

!              Updating max element size
               IF((elbox(2*jj,ii)-elbox(2*jj-1,ii)).gt.maxel(jj))
     2            maxel(jj) = elbox(2*jj,ii)-elbox(2*jj-1,ii)
            ENDDO
      ENDDO
      
!     If the elements are larger, set the box number to the first one where they're bigger
      DO ii =1,nsd
            IF (maxel(ii).gt.sb%step(ii))
     2       n(ii) = diff(ii)/maxel(ii)
      ENDDO
      
!     Here's the final number of sb's in each direction!
      sb%n = n
      print *, n

!     Final size of sb
      sb%step = diff/sb%n

!     Dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%box(sb%n(1)*sb%n(2)*sb%n(3)))

!     These sequences are just for allocating sbdim
      ALLOCATE(seq1(sb%n(3)*sb%n(2)),seq2(sb%n(3)*sb%n(1))
     2   ,seq3(sb%n(2)*sb%n(1)))

      seq1=(/(ii, ii=0, sb%n(2)*sb%n(3)-1, 1)/)*sb%n(1)+1
      cnt2=0
      do ii=1,sb%n(1)*sb%n(3)
            seq2(ii)=ii+cnt2*(sb%n(2)-1)*sb%n(1)
            if (MOD(ii,sb%n(1)).eq.0) cnt2=cnt2+1
      end do
      seq3=(/(ii, ii=0, sb%n(1)*sb%n(2)-1, 1)/)+1

      ! Direction 1
      do ii=1,sb%n(1)
         sb%box(seq1+ii-1)%dim(1) = MINVAL(msh%x(1,:))      !! Check dims here
     2       + sb%step(1)*(ii-1)
      end do

      ! Direction 2
      do ii=1,sb%n(2)
         sb%box(seq2+(ii-1)*sb%n(1))%dim(3) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)
      end do

      ! Direction 3
      do ii=1,sb%n(3)
         sb%box(seq3+(ii-1)*sb%n(1)*sb%n(2))%dim(5)=
     2   MINVAL(msh%x(3,:)) + sb%step(3)*(ii-1)
      end do

      sb%box%dim(2) = sb%box%dim(1) + sb%step(1)
      sb%box%dim(4) = sb%box%dim(3) + sb%step(2)
      sb%box%dim(6) = sb%box%dim(5) + sb%step(3)
      sb%minx(1) = minval(sb%box(:)%dim(1))
      sb%minx(2) = minval(sb%box(:)%dim(3))
      sb%minx(3) = minval(sb%box(:)%dim(5))
      

      ! Finding the sbs the box vertices are in
      DO ii=1,msh%Nel
         elvert(1,1) = elbox(1,ii)
         elvert(2,1) = elbox(3,ii)

         elvert(1,2) = elbox(2,ii)
         elvert(2,2) = elbox(3,ii)

         elvert(1,3) = elbox(2,ii)
         elvert(2,3) = elbox(4,ii)

         elvert(1,4) = elbox(1,ii)
         elvert(2,4) = elbox(4,ii)

         IF (nsd.eq.3) THEN
         
         elvert(3,1) = elbox(5,ii)
         elvert(3,2) = elbox(5,ii)
         elvert(3,3) = elbox(5,ii)
         elvert(3,4) = elbox(5,ii)

         elvert(1,5) = elbox(1,ii)
         elvert(2,5) = elbox(3,ii)
         elvert(3,5) = elbox(6,ii)

         elvert(1,6) = elbox(2,ii)
         elvert(2,6) = elbox(3,ii)
         elvert(3,6) = elbox(6,ii)

         elvert(1,7) = elbox(2,ii)
         elvert(2,7) = elbox(4,ii)
         elvert(3,7) = elbox(6,ii)

         elvert(1,8) = elbox(1,ii)
         elvert(2,8) = elbox(4,ii)
         elvert(3,8) = elbox(6,ii)

         END IF

!        Just doing subroutine idsb, but I haven't made the searchboxes yet so I can't call it
         addloop: DO jj = 1,2**nsd
            ! Set domain back to zero
            xzero(1) = elvert(1,jj) - sb%minx(1)
            xzero(2) = elvert(2,jj) - sb%minx(2)
            xzero(3) = elvert(3,jj) - sb%minx(3)

            ! Find which searchbox the particle is in
            ! Number of searchbox steps in x,y,and z
            xsteps = FLOOR(xzero/sb%step)

            ! Searchbox the element is in
            iSb = xsteps(1) + sb%n(1)*xsteps(2) +
     2     sb%n(1)*sb%n(2)*xsteps(3) + 1

!           Add element to sb (if sb exists)
                  IF ((iSb.gt.0)
     2    .and.(iSb.le.sb%n(1)*sb%n(2)*sb%n(3))) THEN

!           If this isn't the first element going in the box
                  IF (ALLOCATED(sb%box(iSb)%els))THEN
!           First check if the element has been added already
                        IF (ANY(sb%box(iSb)%els.eq.ii))
     2                  cycle addloop

                        ALLOCATE(tmpbox%els(
     2            size(sb%box(iSb)%els)+1))
                        tmpbox%els(1:size(sb%box(iSb)%els))
     2            = sb%box(iSb)%els
                        DEALLOCATE(sb%box(iSb)%els)
                        ALLOCATE(sb%box(iSb)%els(
     2            size(tmpbox%els)))
                        sb%box(iSb)%els = tmpbox%els
                        sb%box(iSb)%els(size(tmpbox%els)) =
     2            ii
                        DEALLOCATE(tmpbox%els)
!           If this is the first element in the box
                  ELSE
                        ALLOCATE(sb%box(iSb)%els(1))
                        sb%box(iSb)%els(1) = ii
                  END IF
                  END IF
         END DO addloop
      ENDDO

!     Same process as above, but got faces
      do ii = 1,msh%nFa
            do jj = 1,msh%fa(ii)%nEl
                  do kk = 1,nsd
                        elboxf(2*kk-1) = 
     2            MINVAL(msh%x(kk,msh%fa(ii)%IEN(:,jj)))
                        elboxf(2*kk  ) =
     2            MAXVAL(msh%x(kk,msh%fa(ii)%IEN(:,jj)))
                  end do

                  elvert(1,1) = elboxf(1)
                  elvert(2,1) = elboxf(3)
         
                  elvert(1,2) = elboxf(2)
                  elvert(2,2) = elboxf(3)
         
                  elvert(1,3) = elboxf(2)
                  elvert(2,3) = elboxf(4)
         
                  elvert(1,4) = elboxf(1)
                  elvert(2,4) = elboxf(4)
         
                  IF (nsd.eq.3) THEN
                  
                  elvert(3,1) = elboxf(5)
                  elvert(3,2) = elboxf(5)
                  elvert(3,3) = elboxf(5)
                  elvert(3,4) = elboxf(5)
         
                  elvert(1,5) = elboxf(1)
                  elvert(2,5) = elboxf(3)
                  elvert(3,5) = elboxf(6)
         
                  elvert(1,6) = elboxf(2)
                  elvert(2,6) = elboxf(3)
                  elvert(3,6) = elboxf(6)
         
                  elvert(1,7) = elboxf(2)
                  elvert(2,7) = elboxf(4)
                  elvert(3,7) = elboxf(6)
         
                  elvert(1,8) = elboxf(1)
                  elvert(2,8) = elboxf(4)
                  elvert(3,8) = elboxf(6)
                  END IF
                  
!        Just doing subroutine idsb, but I haven't made the searchboxes yet so I can't call it
         addloop2: DO kk = 1,2**nsd
            ! Set domain back to zero
            xzero(1) = elvert(1,kk) - sb%minx(1)
            xzero(2) = elvert(2,kk) - sb%minx(2)
            xzero(3) = elvert(3,kk) - sb%minx(3)

            ! Find which searchbox the particle is in
            ! Number of searchbox steps in x,y,and z
            xsteps = FLOOR(xzero/sb%step)

            ! SB el is in
            iSb = xsteps(1) + sb%n(1)*xsteps(2) +
     2     sb%n(1)*sb%n(2)*xsteps(3) + 1

!           Add element to sb (if sb exists/element hasn't been added)
                  IF ((iSb.gt.0)
     2    .and.(iSb.le.sb%n(1)*sb%n(2)*sb%n(3))) THEN
!           First we need to allocate the face structure of the sb...
                        IF(.not.ALLOCATED(sb%box(iSb)%fa))
     2                  ALLOCATE(sb%box(iSb)%fa(msh%nFa))             

!           If this isn't the first element going in the box
                  IF (ALLOCATED(sb%box(iSb)%fa(ii)%els))THEN
!           First check if the element has been added already
                        IF (ANY(sb%box(iSb)%fa(ii)%els.eq.jj))
     2                  cycle addloop2

                        ALLOCATE(tmpbox%els(
     2            size(sb%box(iSb)%fa(ii)%els)+1))
                        tmpbox%els(1:size(sb%box(iSb)%fa(ii)%els))
     2            = sb%box(iSb)%fa(ii)%els
                        DEALLOCATE(sb%box(iSb)%fa(ii)%els)
                        ALLOCATE(sb%box(iSb)%fa(ii)%els(
     2            size(tmpbox%els)))
                        sb%box(iSb)%fa(ii)%els = tmpbox%els
                        sb%box(iSb)%fa(ii)%els(size(tmpbox%els)) =
     2            jj
                        DEALLOCATE(tmpbox%els)
!           If this is the first element in the box
                  ELSE
                        ALLOCATE(sb%box(iSb)%fa(ii)%els(1))
                        sb%box(iSb)%fa(ii)%els = jj
                  END IF
                  END IF
         END DO addloop2
         
            end do
      end do

      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSbe

!-------------------------------------------------------------------- Particles
      SUBROUTINE newSbp(sb,msh,np)
      CLASS(mshType), INTENT(IN) :: msh
      INTEGER, INTENT(IN) :: np
      CLASS(sbpType), INTENT(INOUT):: sb
      INTEGER :: ii,cnt2, iSb(2**nsd),xsteps(2**nsd)
     2  , n(nsd)
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), xzero(nsd), order(nsd,2)
     2 , s
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
      order(2,1) = MINLOC(diff,1, MASK = orderl)
      orderl(order(2,1)) = .false.
      order(3,1) = MINLOC(diff,1, MASK = orderl)
      order(:,2) = diff(order(:,1))
      order(:,2) = order(:,2)/order(1,2)

!     Scaling to get approximately cubic SBs equal to approx number of particles
      s = (np/(order(2,2)*order(3,2)))**(1D0/3D0)
      
!     First n estimate
      DO ii = 1,nsd
            n(order(ii,1)) = INT(s*order(ii,2))
      ENDDO
      
!     Size of sb
      sb%step = diff/((n + 1D0)/2D0)

!     Here's the final number of sb's in each direction, 2 are needed at least !! to be fixed to make sure vel constraint is satisfied
      DO ii = 1,nsd
            sb%n(ii) = MAX(n(ii),2)
      ENDDO

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%box(sb%n(1)*sb%n(2)*sb%n(3)))

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(sb%n(3)*sb%n(2)),seq2(sb%n(3)*sb%n(1))
     2   ,seq3(sb%n(2)*sb%n(1)))

      seq1=(/(ii, ii=0, sb%n(2)*sb%n(3)-1, 1)/)*sb%n(1)+1
      cnt2=0
      do ii=1,sb%n(1)*sb%n(3)
            seq2(ii)=ii+cnt2*(sb%n(2)-1)*sb%n(1)
            if (MOD(ii,sb%n(1)).eq.0) cnt2=cnt2+1
      end do
      seq3=(/(ii, ii=0, sb%n(1)*sb%n(2)-1, 1)/)+1

      ! Allocating sb, such that they overlap by 50%
      ! Direction 1
      do ii=1,sb%n(1)
         sb%box(seq1+ii-1)%dim(1) = MINVAL(msh%x(1,:)) 
     2       + sb%step(1)*(ii-1)/2
      end do

      ! Direction 2
      do ii=1,sb%n(2)
         sb%box(seq2+(ii-1)*sb%n(1))%dim(3) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)/2
      end do

      ! Direction 3
      do ii=1,sb%n(3)
         sb%box(seq3+(ii-1)*sb%n(1)*sb%n(2))%dim(5)=
     2   MINVAL(msh%x(3,:)) + sb%step(3)*(ii-1)/2
      end do

      sb%box%dim(2) = sb%box%dim(1) + sb%step(1)
      sb%box%dim(4) = sb%box%dim(3) + sb%step(2)
      sb%box%dim(6) = sb%box%dim(5) + sb%step(3)
      sb%minx(1) = minval(sb%box(:)%dim(1))
      sb%minx(2) = minval(sb%box(:)%dim(3))
      sb%minx(3) = minval(sb%box(:)%dim(5))

      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSbp
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
      REAL(KIND=8) :: N(msh%eNoN), xl(nsd,msh%eNoN)
      TYPE(pRawType), POINTER :: p
      TYPE(boxelType),  POINTER :: b
      INTEGER :: ii, ind

      N= -1D0
      p => prt%dat(ip)
      IF (p%sbIDe .eq. -1) RETURN
      b => prt%sbe%box(p%sbIDe)

      do ii=1,size(b%els)+1

            IF (ii.eq.1) THEN
                  ind = p%eIDo
            ELSE
                  ind = b%els(ii-1)
            END IF

            xl = msh%x(:,msh%IEN(:,ind))
            N = msh%nAtx(p%x,xl)

            ! Checking if all shape functions are positive
            IF (ALL(N.ge.-1.0D-7)) then
                  p%eID=ind
                  p%N = N
                  EXIT
            END IF
         
      end do

      ! If it loops through everything and doesn't yield a positive shape function,
      ! the particle is outside the domain.

      !! Needs fixing for no tolerancing now (will need to be done in idsbe so it lets us know if we're oob there)
      RETURN
      END FUNCTION shapeFPrt
!--------------------------------------------------------------------
      SUBROUTINE seedPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      INTEGER ip
      TYPE(pRawType) p
      REAL, ALLOCATABLE :: N(:)

      ALLOCATE(N(prt%dmn%msh(1)%eNoN))

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
           !p%x(3) = 0.21D0
           !if (ip.eq.2) p%x(3)=0.1D0
           p%sbIDp = prt%sbp%id(p%x)
           p%sbIDe = prt%sbe%id(p%x)
           prt%dat(ip) = p
           N = prt%shapef(ip,prt%dmn%msh(1))
      END DO

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
      do ii=1,nsd
         do jj=1,msh%eNoN
            fvel(ii) = fvel(ii) + u%v(ii,msh%IEN(jj,p%eID))*p%N(jj)
         end do
      end do

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

      p1 => prt%dat(id1)
      p2 => prt%dat(id2)
      dp  = prt%mat%D()

      p1%OthColl(id2) = .true.
      p2%OthColl(id1) = .true.

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
      IF (ANY(Np1.lt.-1D-7) .or. ANY(Np2.lt.-1D-7)) THEN
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
      prt%collpair(prt%collcnt,:) = (/id1,id2/)
      prt%collt(prt%collcnt) = tcr! + (dt) !something like this                    !! This needs to account for particles that have collided already

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

      ! Note that perpendicular velocities don't change, so we only need to calculate parallel
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
     2 xl(nsd,msh%eNoN)

      p => prt%dat(idp)
      b => prt%sbe%box(p%sbIDe)

      p%faID = 0

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
             END DO

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

            IF (ALL(N.ge.-1D-7).and. (tc.ge.-1D-7)
     2      .and. tc.lt.prt%dt) THEN
                  p%faID(1) = ii
                  p%faID(2) = b%fa(ii)%els(jj)
                  p%ti = tc
                  RETURN
            ENDIF

      ENDDO
      ENDDO faceloop
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
            Ac = prt%dmn%msh(1)%IEN(jj,p%eID)
            prt%twc%v(:,Ac) = prt%twc%v(:,Ac) +
     2      apd*mP/rhoF/prt%wV(Ac)*p%N(a)
      END DO

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
            p%x = msh%x(:,msh%fa(ii)%IEN(1,rndi)) 

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
      END DO

      x = msh%x(:,Ac)      

      DO g = 1, msh%nG
!     First, we want the Jacobian, which (if shpfns are linear) is the same for all gauss points
      IF (g.EQ.1 .OR. .NOT.msh%lShpF) CALL msh%dNdx(g,x,Nx,Jac)
            DO a = 1, msh%eNoN
!           Numerically integrate to get volume of each node
                  effvol(a) = effvol(a) + msh%N(a,g)*Jac*msh%w(g)
            END DO
      END DO

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

!!    For now, I'm just going to do 1st order. It eases a lot of optimization stuff

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
      END DO

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
      IF (ANY(N .le. -1D-7)) THEN
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
      TYPE(mshType), POINTER :: lM
      INTEGER, ALLOCATABLE :: tmpstck(:), N(:)
      REAL(KIND=8):: dtp,maxdtp,sbdt(nsd),dp,taup,rhoP,mu,tim,
     2 P1, P2, rhoF

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

!! Right now: I'm just going to take min of relaxation time and overall, but will need to make it eq%sb so I can only check neighboring eq%sb's

!     Appropriate time step
      dtp = MIN(taup,dt)
!     Subiterations
      subit = 1!FLOOR(dt/dtp)
      dtp = dt/subit
      eq%dt = dtp
      eq%collt = dtp+1

      DO l = 1,subit

      IF (eq%itr .EQ. 0) THEN
!     Reset SB's particles
      DO i = 1,eq%sbp%n(1)*eq%sbp%n(2)*eq%sbp%n(3) !! Still have 1 single sb loop...
            eq%sbp%box(i)%nprt = 0
            IF (ALLOCATED(eq%sbp%box(i)%c))
     2      DEALLOCATE(eq%sbp%box(i)%c)
            ALLOCATE(eq%sbp%box(i)%c(1))
      END DO
      END IF

      DO i=1,eq%n
!     If it's the first iteration, update last velocity, position, info
            IF (eq%itr .EQ. 0) THEN
                  eq%dat(i)%xo = eq%dat(i)%x
                  eq%dat(i)%uo = eq%dat(i)%u
                  eq%dat(i)%eIDo = eq%dat(i)%eID
                  eq%dat(i)%sbIDeo= eq%dat(i)%sbIDe
                  eq%dat(i)%sbIDpo= eq%dat(i)%sbIDp

!           Put particles in SBs
            DO j = 1,2**nsd

            IF ((eq%dat(i)%sbIDp(j).lt. 
     2       eq%sbp%n(1)*eq%sbp%n(2)*eq%sbp%n(3))
     3      .and. (eq%dat(i)%sbIDp(j).gt.0)) then
!     Increase number in box
                  eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt = 
     2            eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt + 1

            IF (.not. ALLOCATED(eq%sbp%box(eq%dat(i)%sbIDp(j))%c))
     2      ALLOCATE(eq%sbp%box(eq%dat(i)%sbIDp(j))%c(1))       

!     Temporarily hold old searchbox values
                  tmpstck = eq%sbp%box(eq%dat(i)%sbIDp(j))%c
                  DEALLOCATE(eq%sbp%box(eq%dat(i)%sbIDp(j))%c)
!     Increase size of container by 1
                  ALLOCATE(eq%sbp%box(eq%dat(i)%sbIDp(j))%
     2             c(eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt))
!     Add to end (a little different if it's the first element)
                  IF (eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt.ne.1) THEN
                  eq%sbp%box(eq%dat(i)%sbIDp(j))%c
     2       (1:eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt-1)=tmpstck
                  eq%sbp%box(eq%dat(i)%sbIDp(j))%c
     2             (eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt) = i
                  ELSE
                        eq%sbp%box(eq%dat(i)%sbIDp(j))%c(1) = i
                  END IF 
            END IF
            ENDDO
            ENDIF  ! for initial iteration

!           Reset other collisions
            IF(.not.ALLOCATED(eq%dat(i)%OthColl))
     2       ALLOCATE(eq%dat(i)%OthColl(eq%n))
            eq%dat(i)%OthColl = .false.

!     Set position and velocity to old variables in preparation for iteration with INS
            eq%dat(i)%x = eq%dat(i)%xo
            eq%dat(i)%u = eq%dat(i)%uo
            eq%dat(i)%eID = eq%dat(i)%eIDo
            eq%dat(i)%sbIDe = eq%dat(i)%sbIDeo
            eq%dat(i)%sbIDp = eq%dat(i)%sbIDpo

!     Set initial advancing time step to solver time step
            eq%dat(i)%remdt = dtp
      ENDDO

      DO i = 1,eq%n
            DO j=1,2**nsd
            IF ((eq%dat(i)%sbIDp(j).gt.0).and.
     2 (eq%dat(i)%sbIDp(j).lt.eq%sbp%n(1)*eq%sbp%n(2)*eq%sbp%n(3))) THEN
                  DO k=1,eq%sbp%box(eq%dat(i)%sbIDp(j))%nprt
                        i2 = eq%sbp%box(eq%dat(i)%sbIDp(j))%c(k)        
!     Check if the particle collides with any other particles during step. Then add to list if so
               IF ((i.lt.i2) .and. .not.(eq%dat(i)%OthColl(i2))) THEN
                        CALL eq%findcoll(i,i2,lM)
               END IF
                  ENDDO
            END IF
            ENDDO
      ENDDO

!     Keep looping over collisions, adding more as they appear, until there aren't any left
      DO WHILE(MINVAL(eq%collt).le.dtp)
!     Enact collisions, starting with shortest time to collision, check for more collisions, add any you find in
            citer = MINLOC(eq%collt, 1)
            i1 = eq%collpair(citer,1)
            i2 = eq%collpair(citer,2)
            CALL eq%collide(i1,i2)
            eq%collt(citer) = dtp + 1
!     Get rid of any previous collisions with these particles
            eq%collt(PACK(eq%collpair
     2 , eq%collpair .eq. i1)) = dtp + 1
            eq%collt(PACK(eq%collpair
     2 , eq%collpair .eq. i2)) = dtp + 1
            eq%dat(i1)%OthColl = .false.
            eq%dat(i2)%OthColl = .false.
            DO i = 1,eq%n                 !! still loops over all other parts if collision occurs
                                          !! to fix, remove these particles from sb, add to new ones
                  IF((i.ne.i1).and.
     2              .not.(eq%dat(i1)%OthColl(i)))
     3 CALL eq%findcoll(i1,i,lM)

                  IF((i.ne.i2).and.
     2              .not.(eq%dat(i2)%OthColl(i)))
     3 CALL eq%findcoll(i2,i,lM)
            ENDDO
      ENDDO

      DO i = 1,eq%n
            CALL eq%adv(i)
            IF (eq%itr .EQ. 0)  print *, eq%dat(i)%x(3),
     2            eq%dat(i)%u(3)
      END DO
            
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

!     Updating norm for solution control
      CALL eq%upNorm(eq%ls%RI%iNorm)
      
!     Checking for exceptions
      CALL io%w%checkException()

      RETURN
      END SUBROUTINE solvePrt
      
      END MODULE PRTMOD


      !! Urgent fixes:
      !! Combine wall into collisions matrix
      !! Collisions still checks every other particle after it has already collided
      !! Mult collisions still need quite a bit of work to keep time consistent (if remdt1 .ne. remdt2)
      !! Don't have it implemented so it can hit multiple walls
