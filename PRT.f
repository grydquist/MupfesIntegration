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
!     Elements contained in searchbox
            INTEGER, ALLOCATABLE :: els(:)
!     Face elements in searchbox (face and element)
            TYPE(facelsType), ALLOCATABLE :: fa(:)
      END TYPE

!     Individual box (particles)
      TYPE boxpType
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
         REAL(KIND=8) :: minx(3)
      CONTAINS
!     Sets up the search boxes pertaining to msh
         PROCEDURE :: new => newSbe
!     Deallocates structure
         PROCEDURE :: free => freeSbe
!     Returns serach box ID, provided the position of a point
         PROCEDURE :: id => idSbe
      END TYPE

!     Search boxes for particle collisions
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
         REAL(KIND=8) :: minx(3)
!     Overlap
         REAL(KIND=8) :: Ov(3)
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
!     Coupling (1,2,4)
         INTEGER :: couple = 1
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
      
!     IMPLEMENTATION

!#####################################################################
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

!     this is just a placeholder for now

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
      CHARACTER(LEN=1), PARAMETER :: names = 'dID'
      INTEGER i, ip, l, fid, tN, pN, n
      INTEGER(KIND=8) pos(2)
      CHARACTER(LEN=stdL) hdr
      TYPE(pRawType), POINTER :: p
      TYPE(prtType) prt, tPrt
      INTEGER, ALLOCATABLE :: lN(:)
      REAL(KIND=8), ALLOCATABLE :: tmpX(:,:), tmpV(:,:), 
     2   tmpI(:)

      tn = s%n
      fid = 142
      OPEN(fid, FILE=TRIM(fName))
      CLOSE(fid,STATUS='DELETE')
      OPEN(fid, FILE=TRIM(fName), ACCESS="STREAM", 
     2      CONVERT='BIG_ENDIAN')
      
      hdr = "# vtk DataFile Version 3.0"//eol//
     1      "Particle info"//eol//
     2      "BINARY"//eol//
     3      "DATASET POLYDATA"//eol//
     4      "POINTS "//STR(tN,intPr)//" double"//eol
      l   = LEN(TRIM(hdr))
      WRITE(fid) hdr(1:l)
!     Four pos point to the position in the file associated with
!     position, velocity, temprature, disID
      pos = l + 1
      pos(2:) = pos(2:) + 24*INT(tN,8) ! position
      hdr = eol//"VERTICES 1 "//STR(tN+1,intPr)//eol
      l   = LEN(TRIM(hdr))
      WRITE(fid,pos=pos(2)) hdr(1:l)
      pos(2:) = pos(2:) + l
      WRITE(fid,pos=pos(2)) tN, (i, i=0,tN-1)
      pos(2:) = pos(2:) + 4*INT(tN,8) + 4 ! Vertices
      hdr = eol//"POINT_DATA "//STR(tN,intPr)//eol
     2      //"VECTORS Vel double"//eol
      l   = LEN(TRIM(hdr))
      WRITE(fid,pos=pos(2)) hdr(1:l)
      pos(2:) = pos(2:) + l

!     Total number of processed particles
      pN = 0
      DO 
            n = tN-pN
            IF (n .EQ. 0) EXIT

            ALLOCATE(tmpV(3,n), tmpX(3,n))
            DO ip=1, n
            i = s%dat(ip)%pid
            tmpX(:,i-pN) = s%dat(ip)%x
            tmpV(:,i-pN) = s%dat(ip)%u
            END DO
            WRITE(fid,pos=pos(1)) tmpX
            WRITE(fid,pos=pos(2)) tmpV
            pos(1:2) = pos(1:2) + 24*n
            DEALLOCATE(tmpV, tmpX)
            pN = pN + n
      END DO
      CLOSE(fid)

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
      INTEGER nFa,iFa, typ2, typ3
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

!     Get the coupling strategy
      typ3 = 0
      lPt1 => lst%get(typ3,"Coupling")
      IF((typ3 .NE. 1) .AND. (typ3 .NE. 2) .AND. 
     2   (typ3 .NE. 4) .AND. (typ3 .NE. 0)) THEN
            io%e = "Simulation must be 1, 2, or 4 way coupled"
      END IF
      IF(typ3 .NE. 0) eq%couple = typ3

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
      CALL eq%sbe%new(eq%dmn%msh(1), eq%n)
      
!     Getting volumes of influence for each node for twc (if tw coupled)
      IF (eq%couple .GT. 1) THEN
            DO i = 1,eq%dmn%msh(1)%nEl
                  volt = effvolel(eq%dmn%msh(1),i)
                  DO a = 1,eq%dmn%msh(1)%eNoN
                        Ac = eq%dmn%msh(1)%IEN(a,i)
                        eq%wV(Ac) = eq%wV(Ac) + volt(a)
                  ENDDO
            ENDDO
      END IF

      CALL eq%seed()

      END SUBROUTINE setupPrt
!---------------------------------------------------------------Elements
      SUBROUTINE newSbe(sb,msh,np)
      CLASS(mshType), INTENT(IN) :: msh
      CLASS(sbeType), INTENT(INOUT):: sb
      INTEGER, INTENT(IN) :: np
      INTEGER :: i, j, k, l, iSb, xsteps(nsd)
     2  , xstepst(nsd,2), iSBmin, SBt, cx, cy, cz, nt
     3  ,  faceplne, elnmfc
      REAL(KIND=8) :: diff(nsd),elvert(nsd,2),xzerot(nsd,2),order(nsd,2)
     2 , s,xzero(nsd), N(msh%eNoN), xl(nsd,msh%eNoN), sbx(nsd)
     3 , sbxmin(nsd), elfcx(nsd,nsd), fnV(nsd), segpt(nsd,2), df(nsd) 
     4 , d1, d2, d3, segV(nsd), ipt(nsd), sbsegs(nsd*2**(nsd-1),nsd,2)
     5 , c1(nsd), c2(nsd), alph, nalph, Calph
      LOGICAL :: orderl(nsd), inbox, fullacc
      INTEGER, ALLOCATABLE :: elfcpts(:,:), elsegs(:,:)

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

!     Scaling to get approximately cubic SBs equal
!     to much greater than # elements
!     s is how much to increasebox size to get approx ratio
!     of # sbs to # el (1/alph)

!     C (Calph) is to be hardcoded in after determination, but about
!     500 is typically pretty good. Nt time steps may need to
!     be estimated. n (nalph) is based off element type as below
      SELECT CASE(msh%eType)
      CASE(eType_TET)
            nalph = 0.28
            Calph = 503
      CASE(eType_BRK)
            nalph = 0.33
            Calph = 50
      CASE DEFAULT
         io%e = "Element type not supported yet"
      END SELECT

      alph = alphind(Calph, nalph, msh%Nel, np, 125)
!     This alph is w.r.t. nonempty boxes, and we need total boxes for s
!     Nb/Nt is approx = Vol_geom/Vol_bounding_box
      alph = alph*msh%vol/(diff(1)*diff(2)*diff(3))*msh%scaleF**nsd

      s = (msh%nEl/alph/(order(2,2)*order(3,2)))**(1D0/3D0)

!     First n estimate
      DO i = 1,nsd
            sb%n(order(i,1)) = MAX(INT(s*order(i,2)),1)
      ENDDO

!     Do we want to get the exact SBs in the elements?
      fullacc = .true.
      
!     Size of sb
      sb%step = diff/sb%n
      nt = sb%n(1)*sb%n(2)*sb%n(3)
      ALLOCATE(sb%box(nt))

      sb%minx(1) = minval(msh%x(1,:))
      sb%minx(2) = minval(msh%x(2,:))
      sb%minx(3) = minval(msh%x(3,:))

!     Making an array of all possible line segments of SB edges
      sbsegs = 0
      IF (nsd.eq.3) THEN
            sbsegs(9:12,1,1) = 1
            sbsegs((/4,5,7,9/),2,1) = 1
            sbsegs((/6,7,8,12/),3,1) = 1
            sbsegs((/2,5,7,8,9,10,11,12/),1,2) = 1
            sbsegs((/1,4,5,6,7,9,10,12/),2,2) = 1
            sbsegs((/3,4,6,7,8,9,11,12/),3,2) = 1
      ELSE
            sbsegs(3:4,1,1) = 1
            sbsegs(2:3,2,1) = 1
            sbsegs(2:3,1,2) = 1
            sbsegs(1:2,2,2) = 1
      END IF

!     Constructing planes for each face of el & making lines from edges
      SELECT CASE(msh%eType)
      CASE(eType_TET)
!           Tet: 4 faces to check
            elnmfc = 4
            ALLOCATE(elfcpts(nsd,elnmfc))
!           Points on face
            elfcpts(:,1) = (/1,2,3/)
            elfcpts(:,2) = (/1,2,4/)
            elfcpts(:,3) = (/1,3,4/)
            elfcpts(:,4) = (/2,3,4/)
!           6 edges to check
            ALLOCATE(elsegs(6,2))
            elsegs(:,1) = (/1,1,1,2,2,3/)
            elsegs(:,2) = (/2,3,4,3,4,4/)

      CASE(eType_BRK)
!           Brick: 6 faces to check
            elnmfc = 6
            ALLOCATE(elfcpts(nsd,elnmfc))
!           Points on face
            elfcpts(:,1) = (/1,2,3/)
            elfcpts(:,2) = (/1,2,5/)
            elfcpts(:,3) = (/1,4,5/)
            elfcpts(:,4) = (/5,6,7/)
            elfcpts(:,5) = (/3,4,7/)
            elfcpts(:,6) = (/2,3,6/)
!           12 edges to check
            ALLOCATE(elsegs(12,2))
            elsegs(:,1) = (/1,1,1,3,3,3,6,6,6,8,8,8/)
            elsegs(:,2) = (/2,4,5,2,4,7,2,5,7,4,5,7/)

      CASE DEFAULT
         io%e = "Element type not supported yet"
      END SELECT

!     Finding the sbs the box vertices are in
      DO i=1,msh%Nel
            xl = msh%x(:,msh%IEN(:,i))
            elvert(1,1) = MINVAL(xl(1,:))
            elvert(1,2) = MAXVAL(xl(1,:))
            elvert(2,1) = MINVAL(xl(2,:))
            elvert(2,2) = MAXVAL(xl(2,:))

         IF (nsd.eq.3) THEN
            elvert(3,1) = MINVAL(xl(3,:))
            elvert(3,2) = MAXVAL(xl(3,:))
         END IF

!           Set domain back to zero
            xzerot(1,:) = elvert(1,:) - sb%minx(1)
            xzerot(2,:) = elvert(2,:) - sb%minx(2)
            xzerot(3,:) = elvert(3,:) - sb%minx(3)

!           Number of sB steps in x,y,and z for both extreme vertices
!           Also, special exception for elvert = exactly max domain val
            xstepst(1,:) = FLOOR(xzerot(1,:)/sb%step(1))
            IF(xstepst(1,1).eq.sb%n(1)) xstepst(1,1) = sb%n(1) - 1
            IF(xstepst(1,2).eq.sb%n(1)) xstepst(1,2) = sb%n(1) - 1

            xstepst(2,:) = FLOOR(xzerot(2,:)/sb%step(2))
            IF(xstepst(2,1).eq.sb%n(2)) xstepst(2,1) = sb%n(2) - 1
            IF(xstepst(2,2).eq.sb%n(2)) xstepst(2,2) = sb%n(2) - 1

            xstepst(3,:) = FLOOR(xzerot(3,:)/sb%step(3))
            IF(xstepst(3,1).eq.sb%n(3)) xstepst(3,1) = sb%n(3) - 1
            IF(xstepst(3,2).eq.sb%n(3)) xstepst(3,2) = sb%n(3) - 1

!           Difference between SB steps for extreme vertices
            xsteps = xstepst(:,2) - xstepst(:,1)

!           Furthest back SB the element is in
            iSBmin = xstepst(1,1) + sb%n(1)*xstepst(2,1)
     2             +  sb%n(1)*sb%n(2)*xstepst(3,1) + 1
            sbxmin = sb%minx + sb%step*xstepst(:,1)

!           Now with this range, we can find all SBs the element is in
!           Total SBs this particle is in
            SBt = (xsteps(1)+1)*(xsteps(2)+1)*(xsteps(3)+1)
            cx = 0
            cy = 0
            cz = 0

!           Loop over all SBs this element is in and add them
            DO j = 1,SBt
!                 First we need to get the current ID of the SB we're in
                  iSB = iSBmin  + cx + cy*sb%n(1)
     2                          + cz*sb%n(1)*sb%n(2)
           
!                 Check if box is in element, first by checking if box
!                 vertices are in the element
!                 If fullacc false, it'l just add all boxes w/o
!                 checking for intersection. If true, it will check.
                  inbox = .true.
                  IF(fullacc) THEN
                  inbox = .false.                  
                  DO k = 1,2**nsd
!                       x sb vertex
                        
                        IF (k.le.4) THEN
                              sbx(1) = sbxmin(1) + sb%step(1)*cx
                        ELSE
                              sbx(1) = sbxmin(1) + sb%step(1)*(cx+1)
                        END IF

!                       y sb vertex
                        IF (ANY(k.eq.(/1,2,5,6/))) THEN
                              sbx(2) = sbxmin(2) + sb%step(2)*cy
                        ELSE
                              sbx(2) = sbxmin(2) + sb%step(2)*(cy+1)
                        END IF

!                       z sb vertex
                        IF ((mod(k,2) .eq. 0.) .and. (nsd.eq.3)) THEN
                              sbx(3) = sbxmin(3) + sb%step(3)*cz
                        ELSEIF (nsd.eq.3) THEN
                              sbx(3) = sbxmin(3) + sb%step(3)*(cz+1)
                        ENDIF

                        N = msh%nAtx(sbx,xl)

!                       Check if vertex is in element
                        IF(ALL(N .ge. -10*EPSILON(N))) THEN
                              inbox = .true.
                              EXIT
                        ENDIF
                  ENDDO

!                 Now check element vertices in box
                  IF(.not. inbox) THEN
                  DO k = 1,msh%eNoN
                       IF((xl(1,k).ge.(sbxmin(1) + sb%step(1)* cx   ))
     2               .and.(xl(1,k).le.(sbxmin(1) + sb%step(1)*(cx+1))) 
     3               .and.(xl(2,k).ge.(sbxmin(2) + sb%step(2)* cy   ))
     4               .and.(xl(2,k).le.(sbxmin(2) + sb%step(2)*(cy+1))))
     5                  THEN
                       IF (nsd.eq.2) THEN 
                        inbox = .true.
                        EXIT
                   ELSEIF((xl(3,k).ge.(sbxmin(3) + sb%step(3)* cz   ))
     2               .and.(xl(3,k).le.(sbxmin(3) + sb%step(3)*(cz+1))))
     3                  THEN
                        inbox = .true.    
                        EXIT
                        ENDIF
                        ENDIF
                  ENDDO
                  ENDIF

!                 Now check if an element edge runs through a box face
                  IF ((nsd.eq.3) .and. .not. inbox) THEN
                  DO k = 1,size(elsegs(:,1))
                        IF (inbox) EXIT
                        segpt(:,1) =  msh%x(:,msh%IEN(elsegs(k,1),i))
                        segpt(:,2) =  msh%x(:,msh%IEN(elsegs(k,2),i))
                        segV = segpt(:,1) - segpt(:,2)
                        DO l =1,2*nsd
                              IF(inbox) EXIT

                              sbx(1) = sbxmin(1) + sb%step(1)*cx
                              sbx(2) = sbxmin(2) + sb%step(2)*cy
                              sbx(3) = sbxmin(3) + sb%step(3)*cz

!                             Faces of SB to check
                              SELECT CASE(l)
                              CASE(1) ! -x
                                    fnV = (/1,0,0/)
                                    faceplne = 1
                              CASE(2) ! +x
                                    fnV = (/1,0,0/)
                                    sbx(1) = sbx(1) +sb%step(1)
                                    faceplne = 1
                              CASE(3) ! -y
                                    fnV = (/0,1,0/)
                                    faceplne = 2
                              CASE(4) ! +y
                                    fnV = (/0,1,0/)
                                    sbx(2) = sbx(2) +sb%step(2)
                                    faceplne = 2
                              CASE(5) ! -z
                                    fnV = (/0,0,1/)
                                    faceplne = 3
                              CASE(6) ! +z
                                    fnV = (/0,0,1/)
                                    sbx(3) = sbx(3) +sb%step(3)
                                    faceplne = 3
                              END SELECT

                              df = segpt(:,1) - sbx
                              d1 = df(1)*fnV(1) + df(2)*fnV(2)
     2                           + df(3)*fnV(3)
                              d2 = segV(1)*fnV(1) + segV(2)*fnV(2)
     2                           + segV(3)*fnV(3)
!                             Parallel if d2 = 0
                              IF (d2 .eq. 0) CYCLE
                              d3 = d1/d2
                              
!                             Intersection pt of line w/ box face plane
                              ipt = segpt(:,1) - segV*d3
!                             Is this point within the segment?
                              IF(d3.lt.0 .or. d3.gt.1) CYCLE
!                             Is this point inside the box face?
                              SELECT CASE(faceplne)
                              CASE(1) ! x face
                                IF((ipt(2).ge.sbx(2))              .and. 
     2                             (ipt(2).le.sbx(2) + sb%step(2)) .and.
     3                             (ipt(3).ge.sbx(3))              .and.
     4                             (ipt(3).le.sbx(3) + sb%step(3))) 
     5                              inbox = .true.
                              CASE(2) ! y face
                                IF((ipt(1).ge.sbx(1))              .and. 
     2                             (ipt(1).le.sbx(1) + sb%step(1)) .and.
     3                             (ipt(3).ge.sbx(3))              .and.
     4                             (ipt(3).le.sbx(3) + sb%step(3))) 
     5                              inbox = .true.
                              CASE(3) ! z face
                                IF((ipt(1).ge.sbx(1))              .and. 
     2                             (ipt(1).le.sbx(1) + sb%step(1)) .and.
     3                             (ipt(2).ge.sbx(2))              .and.
     4                             (ipt(2).le.sbx(2) + sb%step(2))) 
     5                              inbox = .true.
                              END SELECT
                        ENDDO
                  ENDDO
                  ENDIF

!                 SB coords of furthest back vertex
                  sbx(1) = sbxmin(1) + sb%step(1)*cx
                  sbx(2) = sbxmin(2) + sb%step(2)*cy
                  IF (nsd.eq.3) THEN
                        sbx(3) = sbxmin(3) + sb%step(3)*cz
                  ENDIF

!                 Now check if a box edge runs through an element face
                  IF ((nsd.eq.3) .and. .not. inbox) THEN
!                       Loop through all the edges on the SB
                        DO k = 1,nsd*2**(nsd-1)
                              IF(inbox) EXIT
!                             Constructing a line segment w/ faces of SB
                              segpt(1,1) = sbx(1) 
     2                                   + sbsegs(k,1,1)*sb%step(1)
                              segpt(1,2) = sbx(1) 
     2                                   + sbsegs(k,1,2)*sb%step(1)
                              segpt(2,1) = sbx(2) 
     2                                   + sbsegs(k,2,1)*sb%step(2)
                              segpt(2,2) = sbx(2) 
     2                                   + sbsegs(k,2,2)*sb%step(2)
                              segpt(3,1) = sbx(3) 
     2                                   + sbsegs(k,3,1)*sb%step(3)
                              segpt(3,2) = sbx(3) 
     2                                   + sbsegs(k,3,2)*sb%step(3)

!                             Loop through all element faces
                              DO l =1,elnmfc
!                                   Location of face pts
                                    elfcx = msh%x(:,
     2                                      msh%IEN(elfcpts(:,l),i))
!                                   Vector along segment points
                                    segV = segpt(:,1) - segpt(:,2)
!                                   Face normal vector
                                    c1 = elfcx(:,2) - elfcx(:,1)
                                    c2 = elfcx(:,3) - elfcx(:,1)
                                    fnV = cross2(c1,c2)
!                                   Now we can finally find ipt
                                    df = segpt(:,1) - elfcx(:,1)
                                    d1 = df(1)*fnV(1) + df(2)*fnV(2) 
     2                                 + df(3)*fnV(3)
                                    d2 = segV(1)*fnV(1) + segV(2)*fnV(2)
     2                                 + segV(3)*fnV(3)
!                                   Parallel if 0
                                    IF (d2 .eq. 0) CYCLE
                                    d3 = d1/d2
!                                   Intersection pt
                                    ipt = segpt(:,1) - segV*d3
!                                   Within the segment?
                                    IF(d3.lt.0 .or. d3.gt.1) CYCLE
!                                   On element face?
                                    N = msh%nAtx(ipt,xl)
!                                   Perhaps delete?
                                    IF (ALL(N.ge.-1000*EPSILON(N))) THEN
                                          inbox = .true.
                                          EXIT
                                    ENDIF
                              ENDDO
                        ENDDO
                  ENDIF
                  ENDIF

!                 Update counter to go to next SB
                  cx = cx + 1
                  IF (cx .gt. xsteps(1)) THEN
                        cx = 0
                        cy = cy + 1
                  END IF
                  IF (cy .gt. xsteps(2)) THEN
                        cy = 0
                        cz = cz + 1
                  END IF

!                 Add element to sb if sb exists and is in element
                  IF ((iSB.gt.0) .and. (iSB.le.nt) .and. inbox) THEN
!                       If this isn't the first element going in the box
                        IF (ALLOCATED(sb%box(iSB)%els))THEN
                              sb%box(iSB)%els = [sb%box(iSB)%els,i]
!                       If this is the first element in the box
                        ELSE
                              ALLOCATE(sb%box(iSB)%els(1))
                              sb%box(iSB)%els(1) = i
                        END IF
                  END IF
            ENDDO
      ENDDO

!     Faces are just added directly without checking for intersection
      DO i = 1,msh%nFa
            DO j = 1,msh%fa(i)%nEl

                  elvert(1,1) = 
     2             MINVAL(msh%x(1,msh%fa(i)%IEN(:,j)))
                  elvert(1,2) = 
     2             MAXVAL(msh%x(1,msh%fa(i)%IEN(:,j)))
         
                  elvert(2,1) = 
     2             MINVAL(msh%x(2,msh%fa(i)%IEN(:,j)))
                  elvert(2,2) = 
     2             MAXVAL(msh%x(2,msh%fa(i)%IEN(:,j)))
         
                  IF (nsd.eq.3) THEN
                        elvert(3,1) = 
     2             MINVAL(msh%x(3,msh%fa(i)%IEN(:,j)))
                        elvert(3,2) = 
     2             MAXVAL(msh%x(3,msh%fa(i)%IEN(:,j)))
                  END IF
         
!                 Set domain back to zero
                  xzerot(1,:) = elvert(1,:) - sb%minx(1)
                  xzerot(2,:) = elvert(2,:) - sb%minx(2)
                  xzerot(3,:) = elvert(3,:) - sb%minx(3)
         
!                 Find which searchbox the particle is in
!                 Number of searchbox steps in x,y,and z
!                 for both extreme vertices
!                 Also, exception for elvert = exactly max domain value
                  xstepst(1,:) = FLOOR(xzerot(1,:)/sb%step(1))
                  IF(xstepst(1,1).eq.sb%n(1)) xstepst(1,1) = sb%n(1) - 1
                  IF(xstepst(1,2).eq.sb%n(1)) xstepst(1,2) = sb%n(1) - 1

                  xstepst(2,:) = FLOOR(xzerot(2,:)/sb%step(2))
                  IF(xstepst(2,1).eq.sb%n(2)) xstepst(2,1) = sb%n(2) - 1
                  IF(xstepst(2,2).eq.sb%n(2)) xstepst(2,2) = sb%n(2) - 1

                  xstepst(3,:) = FLOOR(xzerot(3,:)/sb%step(3))
                  IF(xstepst(3,1).eq.sb%n(3)) xstepst(3,1) = sb%n(3) - 1
                  IF(xstepst(3,2).eq.sb%n(3)) xstepst(3,2) = sb%n(3) - 1
         
!                 Difference between SB steps for extreme vertices
                  xsteps = xstepst(:,2) - xstepst(:,1)
         
!                 Furthest back SB the element is in
                  iSBmin = xstepst(1,1) + sb%n(1)*xstepst(2,1) +
     2                        sb%n(1)*sb%n(2)*xstepst(3,1) + 1
         
!                 With this range, we can get all SBs the element is in
!                 Total SBs this particle is in
                  SBt = (xsteps(1)+1)*(xsteps(2)+1)*(xsteps(3)+1)
                  cx = 0
                  cy = 0
                  cz = 0

            DO k = 1,SBt
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

!                 Add el to sb (if sb exists/element hasn't been added)
                  IF ((iSB.gt.0).and.(iSB.le.nt)) THEN
!                 First we need to allocate the face structure of the sb
                        IF(.not.ALLOCATED(sb%box(iSB)%fa))
     2                  ALLOCATE(sb%box(iSB)%fa(msh%nFa))             

!                 If this isn't the first element going in the box
                  IF (ALLOCATED(sb%box(iSB)%fa(i)%els))THEN
                        sb%box(iSB)%fa(i)%els =
     2                   [sb%box(iSB)%fa(i)%els,j]
!                 If this is the first element in the box
                  ELSE
                        ALLOCATE(sb%box(iSB)%fa(i)%els(1))
                        sb%box(iSB)%fa(i)%els = j
                  END IF
                  END IF 
            ENDDO       
            ENDDO
      ENDDO

      
! !     Some useful parameters that I used for the paper and outputted

!       d3 = 0 ! mean el/cell (gamma)
!       d2 = 0 ! number of nonzero boxes (Nb, nt is Nbt)
!       DO i = 1,nt
! !           This seems to sort things?
!             do j = size(sb%box(i)%els), 2, -1
!                   call random_number(d1)
!                   cx = int(d1 * j) + 1
!                   cy = sb%box(i)%els(cx)
!                   sb%box(i)%els(cx) = sb%box(i)%els(j)
!                   sb%box(i)%els(j) = cy
!             end do
!             IF (size(sb%box(i)%els) .gt.0) THEN
!                   d3 = d3 + size(sb%box(i)%els)
!                   d2 = d2 + 1
!             END IF
!       ENDDO

!       ! d2 nonzero boxes, d3 (before below) total elements in cells
!       d3 = d3/d2!/nt/3.7412D2*4.8D2
!       ! mean(el/cell), ratio # boxes:els, # boxes, #els, all inc. non0
!       print *, d3, REAL(d2)/REAL(msh%nEl), d2,msh%nEl, nt
!       ! nonsense, alpha, Nbt
!       print *,  REAL(msh%nEl)/REAL(nt), REAL(msh%nEl)/REAL(d2),nt
!       stop

      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSbe

!---------------------------------------------------------------Elements
      SUBROUTINE freeSBe(sb)
      CLASS(sbeType), INTENT(INOUT) :: sb
      
      IF (.not. sb%crtd) RETURN

      DEALLOCATE(sb%box)
      sb%crtd = .false.

      END SUBROUTINE freeSBe

!--------------------------------------------------------------Particles
      SUBROUTINE newSbp(sb,msh,np,dmin)
      CLASS(mshType), INTENT(IN) :: msh
      INTEGER, INTENT(IN) :: np
      REAL(KIND=8), INTENT(IN) :: dmin(nsd)
      CLASS(sbpType), INTENT(INOUT):: sb
      INTEGER :: i
      REAL(KIND=8) :: diff(nsd), step(nsd), s

      IF (sb%crtd) RETURN

!     Domain ranges
      diff(1)=MAXVAL(msh%x(1,:))-MINVAL(msh%x(1,:))
      diff(2)=MAXVAL(msh%x(2,:))-MINVAL(msh%x(2,:))
      diff(3)=MAXVAL(msh%x(3,:))-MINVAL(msh%x(3,:))

!     Starting with two times collision range (min half box)
      step = dmin
      sb%Ov = dmin

!     Approx scaling to get NB = NP
      s = (diff(1)*diff(2)*diff(3)/(np*step(1)*step(2)*step(3)))
     2 **(1D0/3D0)*0.5D0

!     Applying this scaling
      step = step*MAX(s,1D0)
      
!     Getting the number of boxes from this, including min overlap
      sb%n = INT((diff - 2D0*step)/(2D0*step - dmin))
!     Let's get this closer to Np, but keep this ratio
      s = (np/(REAL(sb%n(1)*sb%n(2)*sb%n(3))))**(1D0/3D0)
      sb%n = INT(REAL(sb%n)*s)
      
!     Half-size of sb
      sb%step = (diff + (sb%n-1)*dmin)/(2D0*sb%n)
      sb%nt = sb%n(1)*sb%n(2)*sb%n(3)

      ALLOCATE(sb%box(sb%nt))

      sb%minx(1) = minval(msh%x(1,:))
      sb%minx(2) = minval(msh%x(2,:))
      sb%minx(3) = minval(msh%x(3,:))

!     Allocate SB particles in box
      DO i = 1,sb%nt
            ALLOCATE(sb%box(i)%c(1))
      ENDDO
      sb%box(:)%nprt = 0

      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSbp
!--------------------------------------------------------------Particles
!     Gets rid of SB
      SUBROUTINE freeSBp(sb)
      CLASS(sbpType), INTENT(INOUT):: sb
      
      IF (.not.(sb%crtd)) RETURN

      DEALLOCATE(sb%box)

      sb%crtd = .FALSE.

      END SUBROUTINE freeSBp
!--------------------------------------------------------------Particles
!     Put particles in SBs
      SUBROUTINE addprtSBp(sb,IDSBp,ip)
      CLASS(sbpType), INTENT(INOUT):: sb
      INTEGER, INTENT(IN) :: IDSBp(2**nsd), ip
      INTEGER :: i, cnted(2**nsd)


      cnted = 0
      DO i = 1,2**nsd

!           Actual adding procedure
            IF ((IDSBp(i).le.sb%nt).and. (IDSBp(i).gt.0)
     2       .and. .not.(ANY(IDSBp(i).eq.cnted))) THEN
!                 Sometimes if a particle goers to a new box after collision,
!                 it will keep getting added to the boxes each subiteration.
!                 This is to help prevent that.
                  IF((sb%box(IDSBp(i))%nprt.ne.0) 
     2            .and. ANY(ip .eq. sb%box(IDSBp(i))%c)) CYCLE
!                 Increase number in box
                  sb%box(IDSBp(i))%nprt = sb%box(IDSBp(i))%nprt + 1

!                 Add to end (a little different if it's the first element)
                  IF (sb%box(IDSBp(i))%nprt.ne.1) THEN
                        sb%box(IDSBp(i))%c =
     2                  [sb%box(IDSBp(i))%c,ip] 
                  ELSE
                        DEALLOCATE(sb%box(IDSBp(i))%c)
                        ALLOCATE(sb%box(IDSBp(i))%c(1))
                        sb%box(IDSBp(i))%c = ip
                  END IF
            END IF
            cnted(i) = IDSBp(i)
      ENDDO

      END SUBROUTINE addprtSBp
!--------------------------------------------------------------Particles
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
                  ALLOCATE(sb%box(IDSBp(i))%c(size(tmpstck)))
                  sb%box(IDSBp(i))%c = tmpstck
            END IF
            cnted(i) = IDSBp(i)
      ENDDO

      END SUBROUTINE rmprtSBp
!---------------------------------------------------------------Elements
!     Returns the ID of the searchboxes that contains point x
      FUNCTION idSBe(sb,x) RESULT(iSb)
      IMPLICIT NONE
      CLASS(sbeType), INTENT(IN) :: sb
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER iSb
      REAL(KIND=8) :: xzero(nsd)
      INTEGER :: xsteps(nsd)

!     Set domain back to zero
      xzero(1) = x(1) - sb%minx(1)
      xzero(2) = x(2) - sb%minx(2)
      xzero(3) = x(3) - sb%minx(3)

!     Find which searchbox the particle is in
!     Number of searchbox steps in x,y,and z
      xsteps = FLOOR(xzero/sb%step)
!     Check OOB
      IF (ANY(xsteps.ge.sb%n) .or. ANY(xsteps.lt.0)) THEN
            iSb = -1
            RETURN
      ENDIF
!     SB it's in
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

!     Start all at -1, only positive values are added to boxes
      iSb = 0

!     Set domain back to zero
      xzero = x - sb%minx

!     Find which searchbox the particle is in
!     Number of searchbox steps in x,y,and z, considering overlap
      xsteps =  FLOOR(xzero/(2D0*sb%step-sb%Ov))

!     Set back if it's right on the far edge
      IF(xsteps(1) .eq. sb%n(1)) xsteps(1) = xsteps(1) - 1
      IF(xsteps(2) .eq. sb%n(2)) xsteps(2) = xsteps(2) - 1
      IF(xsteps(3) .eq. sb%n(3)) xsteps(3) = xsteps(3) - 1

!     Furthest searchbox in front
      iSb(1) = xsteps(1) + sb%n(1)*xsteps(2) +
     2   sb%n(1)*sb%n(2)*xsteps(3) + 1

!     Check if still in overlapped portion of prev box, x
      IF(xzero(1) - xsteps(1)*(2D0*sb%step(1) - sb%Ov(1)) 
     2   <= sb%Ov(1) .and. (xsteps(1) .ne. 0)) THEN
            IF (xsteps(1) .ne. 0) iSb(2) = iSb(1) - 1

!           Now check previous sb in y too
            IF(xzero(2) - xsteps(2)*(2D0*sb%step(2) - sb%Ov(2)) 
     2         <= sb%Ov(2) .and. (xsteps(2) .ne. 0)) THEN
                  IF (xsteps(2) .ne. 0) iSb(3) = iSb(1) - sb%n(1)
                  IF (xsteps(1) .ne. 0) iSb(4) = iSb(3) - 1
!                 And z...
                  IF((nsd.eq.3) .and. xzero(3) - 
     2            xsteps(3)*(2D0*sb%step(3) - sb%Ov(3))
     3            <= sb%Ov(3)  .and. (xsteps(3) .ne. 0)) THEN
                        IF (xsteps(3) .ne. 0)
     2                        iSb(5) = iSb(1) - sb%n(1)*sb%n(2)
                        IF (xsteps(1) .ne. 0)
     2                        iSb(6) = iSb(5) - 1
                        IF (xsteps(2) .ne. 0)
     2                        iSb(7) = iSb(5) - sb%n(1)
                        IF (xsteps(1) .ne. 0)
     2                        iSb(8) = iSb(7) - 1
                        RETURN
                  ENDIF
            ENDIF
      ENDIF

!     Not in prev box in x, but still could be in y/z
      IF(xzero(2) - xsteps(2)*(2D0*sb%step(2) - sb%Ov(2)) 
     2     <= sb%Ov(2)) THEN
            IF (xsteps(2) .ne. 0) iSb(3) = iSb(1) - sb%n(1)
!           And z...
            IF((nsd.eq.3) .and. xzero(3) - 
     2       xsteps(3)*(2D0*sb%step(3) 
     3       - sb%Ov(3)) <= sb%Ov(3)) THEN
                  IF (xsteps(3) .ne. 0)
     2                  iSb(5) = iSb(1) - sb%n(1)*sb%n(2)
                  IF (xsteps(2) .ne. 0)
     2                  iSb(7) = iSb(5) - sb%n(1)
                  RETURN
            ENDIF
      ENDIF

!     Now just z
      IF((nsd.eq.3) .and. xzero(3) - 
     2  xsteps(3)*(2D0*sb%step(3) 
     3  - sb%Ov(3)) <= sb%Ov(3)) THEN
            IF (xsteps(3) .ne. 0) iSb(5) = iSb(1) - sb%n(1)*sb%n(2)
            RETURN
      ENDIF

      RETURN
      END FUNCTION idSBp
!--------------------------------------------------------------------
!     Finds the element ID the particle is in and returns shape fns
      FUNCTION shapeFPrt(prt, ip, msh) RESULT(N)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN),TARGET :: prt
      INTEGER, INTENT(IN) :: ip
      TYPE(mshType), INTENT(IN) :: msh
      REAL(KIND=8) :: N(msh%eNoN),xl(nsd,msh%eNoN),Nt(msh%eNoN)
      TYPE(pRawType), POINTER :: p
      TYPE(boxelType),  POINTER :: b
      INTEGER :: i, ind

      N= -1D0
      p => prt%dat(ip)
      IF (p%sbIDe .eq. -1) RETURN
      b => prt%sbe%box(p%sbIDe)
      DO i=1,size(b%els)+1

            IF (i.eq.1) THEN
                  ind = p%eIDo
            ELSE
                  ind = b%els(i-1)
            END IF

            xl = msh%x(:,msh%IEN(:,ind))
            N = msh%nAtx(p%x,xl)

!           Checking if all shape functions are positive
            IF (ALL(N.ge.-EPSILON(N))) then
                  p%eID=ind
                  p%N = N
                  RETURN
            END IF
      ENDDO
         
!     We couldn't find the element. Catch to see if Sb is the issue. Slows
!     things down a lot and should really be run in debug mode.
      DO i = 1,msh%nEl
            xl = msh%x(:,msh%IEN(:,i))
            Nt = msh%nAtx(p%x,xl)

!           Checking if all shape functions are positive
            IF (ALL(Nt.ge.-EPSILON(N))) then
                  print *, i, ip, p%sbIDe,Nt, p%x
                  print *, size(b%els)
                  print *, b%els
                  print *, xl
                  io%e = "SBe issue"
            END IF
      ENDDO

!     If it loops through everything and doesn't yield a positive shape 
!     function, the particle is outside the domain.
      RETURN
      END FUNCTION shapeFPrt
!--------------------------------------------------------------------
!     This could be better. It isn't super random right now.
!     Perhaps could just pick a random location within bounds of domain,
!     and reseed if not in an element? Would be costly, but effective
      SUBROUTINE seedPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      INTEGER ip, iel, i, Npos, temp, cnt
      INTEGER, ALLOCATABLE :: Nind(:)
      TYPE(pRawType) p
      CLASS(mshType), POINTER :: msh
      REAL, ALLOCATABLE :: N(:)
      REAL(KIND=8) randn
      CHARACTER(LEN=stdL) :: fName

      msh => prt%dmn%msh(1)

      ALLOCATE(N(msh%eNoN))
      ALLOCATE(Nind(msh%eNoN))
      Nind = (/1:msh%eNoN/)

      !CALL RSEED(cm%id()) !!
      DO ip=1, prt%n
            p%u    = 0D0
            p%pID  = INT(ip,8)
            CALL RANDOM_NUMBER(randn)
            iel = INT(randn*(msh%nEl - 1),8) + 1
            N = 0
!           Randomizing the order shape functions are put in
            DO i = size(Nind), 2, -1
                  CALL RANDOM_NUMBER(randn)
                  Npos = int(randn*i) + 1
                  temp = Nind(Npos)
                  Nind(Npos) = Nind(i)
                  Nind(i) = temp
            ENDDO
!! Eventually target specific elements using parent domain
            DO i = 1,msh%eNoN - 1
                  cnt = 0
                  DO WHILE ((N(i) .lt. 1D-6) .or. 
     2                    (SUM(N(1:i)) .gt. 1 - (msh%eNoN - i)*1D-6))
                        CALL RANDOM_NUMBER(N(i))
                        N(i) = N(i)*(1 - SUM(N(1:(i - 1))))
                        IF(cnt.eq.99) N = (1/(msh%eNoN+0.001D0))
                        cnt = cnt+1
                  ENDDO
            ENDDO
            N(msh%eNoN) = 1 - SUM(N(1:(msh%eNoN - 1)))
            N = N(Nind)
            N = N/sum(N)

            p%x = 0
            DO i = 1,msh%eNoN
                  p%x = p%x + N(i)*msh%x(:,msh%IEN(i,iel))
            ENDDO

            p%sbIDe = prt%sbe%id(p%x)
!           Reset other collisions
            IF(.not.ALLOCATED(p%OthColl)) ALLOCATE(p%OthColl(0))
            prt%dat(ip) = p
            N = prt%shapef(ip,msh)
      ENDDO
      fName = 
     2    "./alltxtres/prtcsvs/prt_0.vtk"
      CALL writePrt(prt,fName)

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
      INTEGER :: i, j
      TYPE(pRawType), POINTER :: p
!     Derived from flow
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
      DO i=1,nsd
         DO j=1,msh%eNoN
            fvel(i) = fvel(i) + u%v(i,msh%IEN(j,p%eID))*p%N(j)
         ENDDO
      ENDDO

!     Relative velocity
      relvel = fvel - p%u
!     Relative velocity magnitude
      magud = SUM(relvel**2D0)**0.5D0
!     Reynolds Number
      Rep = dp*magud*rhoF/mu
!     Schiller-Neumann (finite Re) correction
      fSN = 1D0 + 0.15D0*Rep**0.687D0
!     Stokes corrected drag force
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
      REAL(KIND=8) :: Np1(nsd+1), Np2(nsd+1), xt1(nsd),xt2(nsd)
      INTEGER, ALLOCATABLE :: tmpcoll(:,:)

      p1 => prt%dat(id1)
      p2 => prt%dat(id2)
!     Weird things happen if both vel are zero (really only at start)
      IF (ALL(abs(p1%u) .lt. EPSILON(p1%u)) .and.
     2    ALL(abs(p2%u) .lt. EPSILON(p2%u))) RETURN
      dp  = prt%mat%D()

      p1%OthColl = [p1%OthColl,id2]
      p2%OthColl = [p2%OthColl,id1]

!     First, check if particles will collide at current trajectory,
!     accounting for if the particle has been advanced for collisions
      a = p1%x(1) - p1%u(1)*(prt%dt - p1%remdt) 
     2  -(p2%x(1) - p2%u(1)*(prt%dt - p2%remdt))
      b = p1%u(1) - p2%u(1)
      c = p1%x(2) - p1%u(2)*(prt%dt - p1%remdt) 
     2  -(p2%x(2) - p2%u(2)*(prt%dt - p2%remdt))
      d = p1%u(2) - p2%u(2)
      if(nsd.eq.3) then
         e = p1%x(3) - p1%u(3)*(prt%dt - p1%remdt) 
     2     -(p2%x(3) - p2%u(3)*(prt%dt - p2%remdt))
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

!     Zeros are when the particle enters or exits vol of another part
      zeros(1) = (-qb + sqrt(qb**2D0-4D0*qa*qc))/(2D0*qa)
      zeros(2) = (-qb - sqrt(qb**2D0-4D0*qa*qc))/(2D0*qa)

!     Negative zeros mean the particle would collide previously in time
      if (ANY(zeros.le.0D0)) RETURN

      tcr = minval(zeros)

!     Exit function if collision won't occur during (remaining)timestep
      IF (tcr.gt.prt%dt) RETURN
!     tcr must also occur after previous collisions with these particles
      IF (tcr.lt.(prt%dt - p1%remdt).or.tcr.lt.(prt%dt - p2%remdt))
     2 RETURN

!     particle locations at point of collision
      xt1 = p1%x
      xt2 = p2%x
!     I change these so I can use shapeF, but I change them back below
      p1%x = p1%u*tcr + p1%x - p1%u*(prt%dt-p1%remdt) 
      p2%x = p2%u*tcr + p2%x - p2%u*(prt%dt-p2%remdt)

!     Check if the particle is outside the domain
      p1%sbIDe = prt%sbe%id(p1%x)
      p2%sbIDe = prt%sbe%id(p2%x)     
      Np1 = prt%shapeF(id1, lM)
      Np2 = prt%shapeF(id2, lM)

!     Change location back
      p1%x = xt1
      p2%x = xt2
      p1%sbIDe = prt%sbe%id(p1%x)
      p2%sbIDe = prt%sbe%id(p2%x)

!     OOB, no collisiion
      IF (ANY(Np1.lt.-EPSILON(Np1)) .or. 
     2    ANY(Np2.lt.-EPSILON(Np2))) THEN
            Np1 = prt%shapeF(id1, lM)
            Np2 = prt%shapeF(id2, lM)
            RETURN
      END IF

      Np1 = prt%shapeF(id1, lM)
      Np2 = prt%shapeF(id2, lM)

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
            prt%collt = [prt%collt,tcr]
!     IF it's not allocated, this is the first collision detected
      ELSE
            ALLOCATE(prt%collt(prt%collcnt))
            prt%collt = tcr
      END IF

      RETURN
      END SUBROUTINE findcollPrt

!--------------------------------------------------------------------
      SUBROUTINE collidePrt(prt,id1,id2,citer)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: citer, id1, id2
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

!     Update location to collision by backtracking then advancing
      p1%x = p1%x - p1%u*(prt%dt - p1%remdt - prt%collt(citer))
      p2%x = p2%x - p2%u*(prt%dt - p2%remdt - prt%collt(citer))

!     Vector parallel and pependicular to collision tangent line
      n1 = (p1%x - p2%x)/((dp + dp)/2)
      n2 = -n1
      temp = cross2(n1,p1%u)
      t1 = cross2(temp,n1)
      temp = cross2(n2,p2%u)
      t2 = cross2(temp,n2)

!     Rare case with no perpendicular velocity
      if (ANY(ISNAN(t1))) t1 = 0D0
      if (ANY(ISNAN(t2))) t2 = 0D0
      
!     Get precollision parallel and perpendicular velocities
      vperp1 = sum(t1*p1%u)
      vpar1  = sum(n1*p1%u)
      vperp2 = sum(t2*p2%u)
      vpar2  = sum(n2*p2%u)

!     Note that perpendicular velocities don't change, 
!     so we only need to calculate parallel
      pa = mp*vpar1 - mp*vpar2
      pb = (-vpar1 - vpar2)*k

      vpar2 = (pa - mp*pb)/(mp + mp)
      vpar1 = pb + vpar2
      vpar2 = -vpar2

!     V here is split into two velocities, so I can add them as vector

      p1%u = vpar1*n1 + vperp1*t1
      p2%u = vpar2*n2 + vperp2*t2

      p1%remdt = prt%dt - prt%collt(citer)
      p2%remdt = prt%dt - prt%collt(citer)

!     Update SB/element/shpfn info info of particles
      ALLOCATE(N(prt%dmn%msh(1)%eNoN))
      IF (prt%couple .EQ. 4) THEN
            p1%sbIDp = prt%sbp%id(p1%x)
            p2%sbIDp = prt%sbp%id(p2%x)
      END IF
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
      CLASS(sbeType), POINTER ::sb
      TYPE(pRawType), POINTER :: p
      TYPE(boxelType),  POINTER :: b
      REAL(KIND=8) :: Jac, xXi(nsd,nsd), Am(nsd,nsd), x1(nsd), tc
      REAL(KIND=8) :: N(msh%eNoN),xi(nsd),Bm(nsd), xc(nsd) 
      INTEGER :: i, j, a,gEl, faceNS(1,msh%fa(1)%eNoN), cnt,k,l
     2 , faceN(msh%fa(1)%eNoN), facev, litr, allb(nsd)
      REAL(KIND=8) s, t, mx, my, ux, uy, uz, lx, ly, lz, iJac,
     2 xl(nsd,msh%eNoN), magud, pV(nsd)

      p  => prt%dat(idp)
      sb => prt%sbe
      b  => sb%box(p%sbIDe)

      p%faID = 0
      magud = SUM(p%u**2D0)**0.5D0

!     First, we need to find all element boxes the particle travels thru
!     pV is the vector detailing the distance the particle travels from
!     the lower corner of the SB it is currently in
      pV = p%u*p%remdt + MOD(p%x,sb%step)
      
!     Find how many SB the particle travels in each direction !!not done
      allb = INT(pV/sb%step)
      faceloop: DO i=1,msh%nFa
      litr = 0
      j = 0
      DO WHILE(litr.eq.0)
            IF (ALLOCATED(b%fa) .and. ALLOCATED(b%fa(i)%els)) THEN !here
                  j = j + 1
            ELSE
                  EXIT
            END IF

            IF (j .eq. size(b%fa(i)%els)) litr = 1
!           First we need to find the volumetric element associated
!           with the face, and get x's associated with that
            gEl = msh%fa(i)%gE(b%fa(i)%els(j))
            xl = msh%x(:,msh%IEN(:,gEl))
!           Next, we find which face of the volumetric element is the 
!           face element
            out: DO  k= 1, msh%fa(i)%eNoN
                  cnt = 1
                  DO l = 1,msh%eNoN 
                        IF (msh%IEN(l,gEl) .eq.
     2                      msh%fa(i)%IEN(k,b%fa(i)%els(j))) THEN
                              faceNS(1,k) = cnt
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
                  xXi(:,1) = xXi(:,1) + xl(:,a) * msh%Nx(1,a,1)
                  xXi(:,2) = xXi(:,2) + xl(:,a) * msh%Nx(2,a,1)
                  xXi(:,3) = xXi(:,3) + xl(:,a) * msh%Nx(3,a,1)   
!                 Location of Gauss point (for Bm)
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

!     Now, knowing the face, we have information about the value of one
!     parent coordinate, which we can use to solve for tc. Then we can
!     find the particle at point of collision with the plane in which
!     the wall face is in. This is the same concept as NAtx, except
!     we're imposing the condition that xc will be on the same plane as
!     the face we're checking

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

!          +y
            CASE(1)
      xi(2) = 1D0
      tc = -(-xi(2)+Am(2,1)*p%x(1) + Am(2,2)*p%x(2) + Am(2,3)*p%x(3)
     2   +  Bm(2))/(Am(2,1)*p%u(1) + Am(2,2)*p%u(2) + Am(2,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)

!           +z
            CASE(2)
      xi(3) = 1D0
      tc = -(-xi(3)+Am(3,1)*p%x(1) + Am(3,2)*p%x(2) + Am(3,3)*p%x(3)
     2   +  Bm(3))/(Am(3,1)*p%u(1) + Am(3,2)*p%u(2) + Am(3,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

!           +x
            CASE(3)
      xi(1) = 1D0
      tc = -(-xi(1)+Am(1,1)*p%x(1) + Am(1,2)*p%x(2) + Am(1,3)*p%x(3)
     2   +  Bm(1))/(Am(1,1)*p%u(1) + Am(1,2)*p%u(2) + Am(1,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

!           -y
            CASE(4)
      xi(2) = -1D0
      tc = -(-xi(2)+Am(2,1)*p%x(1) + Am(2,2)*p%x(2) + Am(2,3)*p%x(3)
     2   +  Bm(2))/(Am(2,1)*p%u(1) + Am(2,2)*p%u(2) + Am(2,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)

!           -z
            CASE(5)
      xi(3) = -1D0
      tc = -(-xi(3)+Am(3,1)*p%x(1) + Am(3,2)*p%x(2) + Am(3,3)*p%x(3)
     2   +  Bm(3))/(Am(3,1)*p%u(1) + Am(3,2)*p%u(2) + Am(3,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)

!           -x
            CASE(6)
      xi(1) = -1D0
      tc = -(-xi(1)+Am(1,1)*p%x(1) + Am(1,2)*p%x(2) + Am(1,3)*p%x(3)
     2   +  Bm(1))/(Am(1,1)*p%u(1) + Am(1,2)*p%u(2) + Am(1,3)*p%u(3))
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
            tc = (1D0-p%x(1) * (Am(1,1)+Am(2,1)+Am(3,1))-
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
!     Having some issues with machine precision. This is to fix that
      xi(1) = xi(1) + (1D0 - sum(xi))
          
            CASE(2)
            xi(3) = 0D0
      tc = -(-xi(3)+Am(3,1)*p%x(1) + Am(3,2)*p%x(2) + Am(3,3)*p%x(3)
     2   +  Bm(3))/(Am(3,1)*p%u(1) + Am(3,2)*p%u(2) + Am(3,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)    

            CASE(3)
            xi(2) = 0D0
      tc = -(-xi(2)+Am(2,1)*p%x(1) + Am(2,2)*p%x(2) + Am(2,3)*p%x(3)
     2   +  Bm(2))/(Am(2,1)*p%u(1) + Am(2,2)*p%u(2) + Am(2,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(1) = Am(1,1)*xc(1) + Am(1,2)*xc(2) + Am(1,3)*xc(3) + Bm(1)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)

            CASE(4)
            xi(1) = 0D0
      tc = -(-xi(1)+Am(1,1)*p%x(1) + Am(1,2)*p%x(2) + Am(1,3)*p%x(3)
     2   +  Bm(1))/(Am(1,1)*p%u(1) + Am(1,2)*p%u(2) + Am(1,3)*p%u(3))
      xc = (p%x + p%u*tc)
      xi(3) = Am(3,1)*xc(1) + Am(3,2)*xc(2) + Am(3,3)*xc(3) + Bm(3)
      xi(2) = Am(2,1)*xc(1) + Am(2,2)*xc(2) + Am(2,3)*xc(3) + Bm(2)  

      END SELECT

            N(1) = xi(1)
            N(2) = xi(2)
            N(3) = xi(3)
            N(4) = 1D0 - xi(1) - xi(2) - xi(3)
!           Machine precision issues still
            IF (facev .eq. 1) N(4) = 0D0
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

!     End NatxiEle

!           Also some parts are getting checked here multiple times per iter??
            IF (ALL(N.ge.-EPSILON(N)) .and. (tc.ge.-EPSILON(N)) 
     2                                .and. (tc.le.p%remdt)) THEN
                  p%faID(1) = i
                  p%faID(2) = b%fa(i)%els(j)
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
      INTEGER :: i, j, rndi, Ac
      REAL(KIND=8) :: dp, rho, k, nV(nsd), tV(nsd),
     2 a(nsd), b(nsd), vpar, vperp, temp(nsd), rnd,
     3 apd(nsd), mp, rhoF

      p   =>prt%dat(idp)
      dp  = prt%mat%D()
      rho = prt%mat%rho()
      k   = prt%mat%krest()
      rhoF  = prt%mns%rho()
      mp = pi*rho/6D0*dp**3D0

!     Advance to collision location
      p%x  = p%u*p%ti + p%x
      p%remdt = p%remdt - p%ti

      IF (faTyp(p%faID(1)) .EQ. 3) THEN
!     Hitting wall

!           Get normal/tangent vector
            a = msh%x(:,msh%fa(p%faID(1))%IEN(1,p%faID(2))) - 
     2          msh%x(:,msh%fa(p%faID(1))%IEN(2,p%faID(2)))
            b = msh%x(:,msh%fa(p%faID(1))%IEN(1,p%faID(2))) - 
     2          msh%x(:,msh%fa(p%faID(1))%IEN(3,p%faID(2)))
            nV = CROSS2(a,b)
            temp = CROSS2(nV,p%u)
            tV = CROSS2(temp,nV)
!          Rare case with no perpendicular velocity
            IF (ANY(ISNAN(tV))) tV = 0D0
            vperp = sum(tV*p%u)
            vpar  = sum(nV*p%u)*k
!          Change velocities to after collision
            p%u = -vpar*nV +vperp*tV

      ELSE

!     Exiting domain
      faloop:DO i = 1,msh%nFa
!           Search for inlet to put particle back into
            IF (faTyp(i) .EQ. 1) THEN
!           Select random node on face to set as particle position
            CALL RANDOM_NUMBER(rnd)
            rndi = FLOOR(msh%fa(i)%nEl*rnd + 1)
            p%x = 0
            DO j = 1, msh%fa(i)%eNoN
                  p%x = p%x + msh%x(:,msh%fa(i)%IEN(j,rndi))
            END DO
!           Put right in middle of face element
            p%x = p%x/msh%fa(i)%eNoN
            EXIT faloop
            END IF
      ENDDO faloop

      END IF

      END SUBROUTINE wallPrt

!--------------------------------------------------------------------
!     Finds the volume influenced by moving one node in one element,
!     for all nodes in that element
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
!     First, we want the Jacobian, which (if shpfns are linear) is the 
!     same for all gauss points
      IF (g.EQ.1 .OR. .NOT.msh%lShpF) CALL msh%dNdx(g,x,Nx,Jac)
            DO a = 1, msh%eNoN
!                 Numerically integrate to get volume of each node
                  effvol(a) = effvol(a) + msh%N(a,g)*Jac*msh%w(g)
            ENDDO
      ENDDO

      END FUNCTION

!--------------------------------------------------------------------
!     I use cross products a couple times above
!     Also normalizes to unit vector
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
      TYPE(mshType), POINTER :: msh
      TYPE(pRawType), POINTER :: p
      TYPE(sbpType), POINTER :: sbp
      TYPE(sbeType), POINTER :: sbe
      REAL(KIND=8) rhoF, rhoP,
     2   mu, mp, g(nsd), dp
      REAL(KIND=8), ALLOCATABLE :: N(:), Ntmp(:)
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd),
     2  tmpwr(nsd),mom, xt(nsd), ut(nsd)
      INTEGER :: a, Ac, i, tsbIDp(2**nsd),tsbIDe, teID

!     Gravity
      g(1) = prt%mat%fx()
      g(2) = prt%mat%fy()
      IF (nsd .eq. 3) g(3) = prt%mat%fz()
      msh => prt%dmn%msh(1)
      ALLOCATE(N(msh%eNoN), Ntmp(msh%eNoN))
      p  => prt%dat(idp)
      IF (prt%couple .EQ. 4) sbp => prt%sbp
      sbe => prt%sbe
!     Original vel and position
      xt = p%x
      ut = p%u
      rhoF  = prt%mns%rho()
      rhoP  = prt%mat%rho()
      mu    = prt%mns%mu()
      dp    = prt%mat%D()
      mp = pi*rhoP/6D0*dp**3D0

!     Get drag acceleration
      apd = prt%drag(idp)
!     Total acceleration (just drag and buoyancy now)
      apT = apd + g*(1D0 - rhoF/rhoP)

!     2nd order advance (Heun's Method)
!     Predictor
      p%x = p%x + p%remdt*p%u
      p%u = p%u + p%remdt*apT

!     Find which searchbox prediction is in
      tsbIDe  = p%sbIDe
      p%sbIDe = sbe%id(p%x)
      N = p%N
!     Get shape functions/element of prediction
      Ntmp = prt%shapeF(idp, msh)
!     Check if predictor OOB. Do first order if so
      IF (ANY(Ntmp.le.-EPSILON(N))) THEN
!           Resets position back, does first order to/beyond wall
            apTpred = apT
      ELSE
!           Get drag acceleration of predicted particle
            apdpred = prt%drag(idp)
            apTpred = apdpred + g*(1D0 - rhoF/rhoP)
      END IF
      p%N = N
      p%sbIDe = tsbIDe

!     Corrector
      p%u = ut + p%remdt*(apT + apTpred)*0.5D0
      p%x = xt + p%remdt*p%u

!     Send drag to fluid (if tw coupled)
      IF(prt%couple .GT. 1) THEN 
            DO a=1,msh%eNoN
                  Ac = msh%IEN(a,p%eID)
                  prt%twc%v(:,Ac) = prt%twc%v(:,Ac) +
     2            0.5D0*(apd + apdpred)*mP/rhoF/prt%wV(Ac)*p%N(a)
            ENDDO
      END IF

!     For checking if momentum in by particles = momentum gain in system
      IF (prt%itr .EQ. 0) THEN
            !tmpwr = apd
            !write(88,*) sqrt(tmpwr(1)**2+tmpwr(2)**2+tmpwr(3)**2)*mp
            !print *, sqrt(tmpwr(1)**2+tmpwr(2)**2+tmpwr(3)**2)*mp
            !mom =prt%dmn%msh(1)%integ(u%v, 3)
            !print *, mom
            !call sleep(1)
      END IF

!     Check if particle went out of bounds continuously
      N = -1D0
      DO WHILE(ANY(N .le. EPSILON(N)))
            tsbIDe = p%sbIDe
            tsbIDp = p%sbIDp
            teID = p%eID
            IF (prt%couple .EQ. 4) p%sbIDp = sbp%id(p%x)
            p%sbIDe = sbe%id(p%x)
            N = prt%shapeF(idp, msh)

!           If OOB, do a wall collision and return to top
            IF (ANY(N .le. -EPSILON(N))) THEN
                  p%x = xt
                  p%sbIDe = tsbIDe
                  p%sbIDp = tsbIDp
                  p%eID = teID
                  CALL prt%findwl(idp,msh)
                  CALL prt%wall(idp,msh)
                  xt = p%x
                  IF (prt%couple .EQ. 4) p%sbIDp = sbp%id(p%x)
                  p%sbIDe = sbe%id(p%x)
                  N = prt%shapeF(idp, msh)
                  N = -1
                  p%x = p%x + p%remdt*p%u
            END IF
      ENDDO

      RETURN

      END SUBROUTINE advPrt
!--------------------------------------------------------------------
      SUBROUTINE solvePrt(eq)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT):: eq
      INTEGER ip, a, i ,j, k,l, subit, citer,i2,i1
     2 , IDSBp(2**nsd), IDSBp1(2**nsd), IDSBp2(2**nsd), i12, i22
      TYPE(mshType), POINTER :: lM
      INTEGER, ALLOCATABLE :: N(:) !!  Get rid of me!!
      REAL(KIND=8):: dtp, dp, taup, rhoP, mu,
     2 P1, P2, rhoF, umax(3) = 0D0, dmin(nsd)
      CHARACTER (LEN=stdl) fName

      lM => eq%dmn%msh(1)
      rhoP  = eq%mat%rho()
      rhoF  = eq%mns%rho()
      mu    = eq%mns%mu()
      dp    = eq%mat%D()
      ALLOCATE(N(lm%eNoN))

!     Reset twc force to zero
      eq%twc%v(:,:) = 0D0

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

!     Collision search boxes: only needed if 4 way coupled
      IF(eq%couple .EQ. 4) THEN
      IF (eq%itr .EQ. 0) THEN
!     Reset SB for particles
            DO i = 1, eq%n
                  umax = MAX(umax,ABS(eq%dat(i)%u))
            ENDDO
            dmin = 2D0*umax*dtp+dp
            CALL eq%sbp%free()
            CALL eq%sbp%new(eq%dmn%msh(1), eq%n,dmin)
      END IF
      END IF

      DO i=1,eq%n
!           If the first iteration, update last velocity, position, info
            IF (eq%itr .EQ. 0) THEN
                  eq%dat(i)%xo = eq%dat(i)%x
                  eq%dat(i)%uo = eq%dat(i)%u
                  eq%dat(i)%eIDo = eq%dat(i)%eID
                  eq%dat(i)%sbIDeo= eq%dat(i)%sbIDe

!                 Put particles in SBs (4-way coupled)
                  IF (eq%couple .EQ. 4) THEN
                        eq%dat(i)%sbIDpo = eq%dat(i)%sbIDp
                        eq%dat(i)%sbIDp = eq%sbp%id(eq%dat(i)%x)
                        idSBp = eq%dat(i)%sbIDp
                        CALL eq%sbp%addprt(idSBp,i)
                  END IF
            ENDIF

!           Set position and velocity to old variables in preparation 
!           for iteration with INS
            eq%dat(i)%x = eq%dat(i)%xo
            eq%dat(i)%u = eq%dat(i)%uo
            eq%dat(i)%eID = eq%dat(i)%eIDo
            eq%dat(i)%sbIDe = eq%dat(i)%sbIDeo

!           Set initial advancing time step to solver time step
            eq%dat(i)%remdt = dtp
!           Reset collisions from last time step
            IF (eq%couple .EQ. 4) THEN
                  DEALLOCATE(eq%dat(i)%OthColl)
                  ALLOCATE  (eq%dat(i)%OthColl(0))
                  eq%dat(i)%sbIDp = eq%dat(i)%sbIDpo
            ENDIF
      ENDDO

!     Collisions: only if 4-way coupled
      IF( eq%couple .EQ. 4) THEN
!     Go through each box and check collisions with all prts in box
      DO i = 1,eq%sbp%nt
            DO j=1,eq%sbp%box(i)%nprt
                  DO k=j+1,eq%sbp%box(i)%nprt
                        i1 = eq%sbp%box(i)%c(j)
                        i2 = eq%sbp%box(i)%c(k)        
!                       Check if the particle collides with any other 
!                       particles during step. Add to list if so
                        IF (.not.ANY(eq%dat(i1)%OthColl.eq.i2)) THEN
                              CALL eq%findcoll(i1,i2,lM)
                        ENDIF
                  ENDDO
            ENDDO
      ENDDO

      IF (ALLOCATED(eq%collt)) THEN
!     Keep looping over collisions, adding more as they appear, 
!     until there aren't any left
      DO WHILE(MINVAL(eq%collt).le.dtp)
!           Enact collisions, starting with shortest time to collision, 
!           check for more collisions, add any you find in
            citer = MINLOC(eq%collt, 1)
            i1 = eq%collpair(citer,1)
            i2 = eq%collpair(citer,2)

!           Remove particles from sb before collision
            idSBp = eq%dat(i1)%sbIDpo
            CALL eq%sbp%rmprt(idSBp,i1)
            idSBp = eq%dat(i2)%sbIDpo
            CALL eq%sbp%rmprt(idSBp,i2)

!           Enact collisions
            CALL eq%collide(i1,i2,citer)

!           Put particles into their new SBs
            idSBp = eq%dat(i1)%sbIDp
            CALL eq%sbp%addprt(idSBp,i1)
            idSBp = eq%dat(i2)%sbIDp
            CALL eq%sbp%addprt(idSBp,i2)

            eq%collt(citer) = dtp + 1
!           Get rid of any previous collisions with these particles
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
!                       Check if the particle collides with other parts
!                       after initial collision. Add to list if so
                        IF ((i1.ne.i12) .and. 
     2                   .not. ANY(eq%dat(i1)%OthColl.eq.i12))
     3                   CALL eq%findcoll(i1,i12,lM)
                  ENDDO
            END IF

            IF ((IDSBp2(i).gt.0).and. (IDSBp2(i).le.eq%sbp%nt)) THEN
                  DO k=1,eq%sbp%box(IDSBp2(i))%nprt
                        i22 = eq%sbp%box(IDSBp2(i))%c(k)         
!                       Check if the particle collides with other parts
!                       after initial collision. Add to list if so
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
!     End of collisions
      END IF

      DO i = 1,eq%n
!           Now no particles will collide on their path. Safe to advance
            CALL eq%adv(i)
      ENDDO
            
      P1 = eq%dmn%msh(1)%integ(1,eq%Pns%s)
      P2 = eq%dmn%msh(1)%integ(2,eq%Pns%s)

      IF (eq%itr .EQ. 0) write(88,*) eq%dat(1)%u(3), !eq%dat(2)%u
     2 eq%dmn%msh(1)%integ(eq%Uns%v, 3)*rhoF,
     3 time
!     3  (eq%dmn%msh(1)%integ(1, eq%Uns%v,3 ) -
!     4  eq%dmn%msh(1)%integ(2, eq%Uns%v,3 ) -
!     5  eq%dmn%msh(1)%integ(3, eq%Uns%v,3 ))*1.2D0,
!     6  ,P1,P2
      ENDDO

!     Write particle data if it's a multiple of ten timestep
      IF (mod(cTs,10).eq.0 .and. eq%itr.eq.0) THEN
            fName = 
     2       "./alltxtres/prtcsvs/prt_"//
     3       STR(cTs)//".vtk"
            CALL writePrt(eq,fName)
      ENDIF

!     Updating norm for solution control
      CALL eq%upNorm(eq%ls%RI%iNorm)
      
!     Checking for exceptions
      CALL io%w%checkException()

      RETURN
      END SUBROUTINE solvePrt
      !--------------------------------------------------------------------
      FUNCTION alphind(C, n, Nelt, Npt, Ntt) RESULT(alph)
!     Secant method to  find alpha opt given C, Nt, Nel, Np, n
      INTEGER :: Nelt, Npt, Ntt
      REAL(KIND=8) :: C, n, Nel, Np, Nt, alph,
     2      x1, x2, f1, f2, Cf, xn, fn
      Nel = REAL(Nelt)
      Np = REAL(Npt)
      Nt = REAL(Ntt)

      Cf = C*Nel/Np/Nt
      x1 = 0
      x2 = 100
      f1 = Cf - x1**(n+1)*(x1**n+1)**(1/n-1)
      f2 = Cf - x2**(n+1)*(x2**n+1)**(1/n-1)
      fn = 1
      
      DO WHILE (abs(fn) .gt. 1e-6)
            xn = x1 - f1/(f1 - f2)*(x1-x2)
            fn = Cf - xn**(n+1)*(xn**n+1)**(1/n-1)
            x2 = x1
            x1 = xn
            f2 = f1
            f1 = fn
      ENDDO
      alph = xn
      END FUNCTION alphind

      
      END MODULE PRTMOD


      !! Urgent fixes:
      !! Doesn't check collisions after wall (maybe make wall collisions part of coll list?)
      !! Wall collisions just use first wall it hits
            ! Very often not a problem, but could be for some domains, including cylindrical
      !! Shapefprt should really be a subroutine, and N should be a part of dat
      !! When a prt hits a wall, check all sbe it passes thru (started