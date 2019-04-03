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

!     Face element boxes
      TYPE facelsboxType
            REAL(KIND=8), ALLOCATABLE :: elbox(:,:)
      END TYPE

!     Faces and their elements
      TYPE facelsType
      ! Elements of face
            INTEGER, ALLOCATABLE :: els(:)
      END TYPE

!     Individual box
      TYPE boxType
!     Box dimensions (minx,maxx,miny,maxy,minz,maxz)
            REAL(KIND=8) :: dim(6)
!     Elements contained in searchbox
            INTEGER, ALLOCATABLE :: els(:)
!     Face elements in searchbox (face and element)
            TYPE(facelsType), ALLOCATABLE :: fa(:)
      END TYPE

!     Search box collection type
      TYPE sbType
!     Created?
         LOGICAL :: crtd = .FALSE.
!     Number of boxes in each direction
         INTEGER n(3)
!     Size of boxes
         REAL(KIND=8) :: step(3)
!     Individual boxes
         TYPE(boxtype), ALLOCATABLE :: box(:)
      CONTAINS
!     Sets up the search boxes pertaining to msh
         PROCEDURE :: new => newSb
!     Returns serach box ID, provided the position of a point
         PROCEDURE :: id => idSb
      END TYPE

!     Your data type that you need to store in a Dynamic Sized Container
      TYPE pRawType
!     Eulerian element ID that this particle belongs to
         INTEGER(KIND=8) :: eID=1
!     Previous element ID
         INTEGER(KIND=8) :: eIDo=1
!     Particle ID
         INTEGER(KIND=8) :: pID=0
!     Searchbox ID
         INTEGER(KIND=8) :: sbID(8) = 0
!     Previous sbID
         INTEGER(KIND=8) :: sbIDo(8) = 0
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
!     Has this particle collided during the current time step?
         LOGICAL :: collided = .FALSE.
!     Remaining time in timestep after collision
         REAL(KIND=8) :: remdt
!     Time until collision
         REAL(KIND=8) :: ti
!     Did particle hit a wall?
         LOGICAL :: wall = .TRUE.
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
         TYPE(sbType) :: sb
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

      eq%wV = 0D0
      cnt =0
      DO i=1, eq%n
         eq%ptr(i) = i
         ALLOCATE(eq%dat(i)%N(eq%dmn%msh(1)%eNoN))
      END DO
      eq%mat  => FIND_MAT('Particle')
      eq%var(1) = gVarType(nsd,'PVelocity',eq%dmn)
      eq%var(2) = gVarType(1,'PPosition',eq%dmn)

!     Assuming everything is in first mesh for searchboxes !!!!!!
      CALL eq%sb%new(eq%dmn%msh(1))
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
!--------------------------------------------------------------------
      SUBROUTINE newSb(sb,msh)
      CLASS(mshType), INTENT(IN) :: msh
      CLASS(sbType), INTENT(INOUT):: sb
      TYPE(facelsboxType), ALLOCATABLE :: facels(:)
      INTEGER :: ii,jj,cnt2,kk,ll
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), elbox(2*nsd,msh%nEl), eps(nsd)
      INTEGER, ALLOCATABLE :: sbel(:),sbelf(:)
      INTEGER :: split(nsd)
      IF (sb%crtd) RETURN

 !!     ! Will need to replace eventually with something that targets element size

      split=(/3,3,3/)
      sb%n=split

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%box(sb%n(1)*sb%n(2)*sb%n(3)))
      ALLOCATE(sbel(msh%nEl))
      ALLOCATE(facels(msh%nFa))
      ALLOCATE(sbelf(MAXVAL(msh%fa(:)%nEl)))

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(sb%n(3)*sb%n(2)),seq2(sb%n(3)*sb%n(1))
     2   ,seq3(sb%n(2)*sb%n(1)))

      ! Domain ranges
      diff(1)=MAXVAL(msh%x(1,:))-MINVAL(msh%x(1,:))
      diff(2)=MAXVAL(msh%x(2,:))-MINVAL(msh%x(2,:))
      diff(3)=MAXVAL(msh%x(3,:))-MINVAL(msh%x(3,:))

      ! Size of sb
      sb%step=(diff)/((sb%n+1D0)/2D0) 
            
      ! Tolerance, based off max between searchbox size (for particle going out of bounds)
      ! and max size needed for id'ing algorithm in idSB
      eps = MAXVAL((/2*diff/(sb%n-1),sb%step/))
      sb%step=(diff+eps)/((sb%n+1D0)/2D0) 
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
     2       + sb%step(1)*(ii-1)/2 - eps(1)/2
      end do

      ! Direction 2
      do ii=1,sb%n(2)
         sb%box(seq2+(ii-1)*sb%n(1))%dim(3) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)/2 - eps(2)/2
      end do

      ! Direction 3
      do ii=1,sb%n(3)
         sb%box(seq3+(ii-1)*sb%n(1)*sb%n(2))%dim(5)=
     2   MINVAL(msh%x(3,:)) + sb%step(3)*(ii-1)/2 - eps(3)/2
      end do

      sb%box%dim(2) = sb%box%dim(1) + sb%step(1)
      sb%box%dim(4) = sb%box%dim(3) + sb%step(2)
      sb%box%dim(6) = sb%box%dim(5) + sb%step(3)

      ! Making boxes surrounding elements
      do ii=1,msh%Nel
         do jj=1,nsd
            elbox(2*jj-1,ii) = MINVAL(msh%x(jj,msh%IEN(:,ii)))
            elbox(2*jj  ,ii) = MAXVAL(msh%x(jj,msh%IEN(:,ii)))
         end do
      end do

!     Make boxes around FACE elements
      do ii = 1,msh%nFa
            ALLOCATE(facels(ii)%elbox(2*nsd,msh%fa(ii)%nEl))
            do jj = 1,msh%fa(ii)%nEl
                  do kk = 1,nsd
                        facels(ii)%elbox(2*kk-1,jj) = 
     2            MINVAL(msh%x(kk,msh%fa(ii)%IEN(:,jj)))
                        facels(ii)%elbox(2*kk  ,jj) =
     2            MAXVAL(msh%x(kk,msh%fa(ii)%IEN(:,jj)))
                  end do
            end do
      end do

!     Populate boxes with elements
      do ii = 1,sb%n(1)*sb%n(2)*sb%n(3)
         ALLOCATE(sb%box(ii)%fa(msh%nFa))
         cnt2 = 1
         sbel = 0
         outer: do jj = 1,msh%Nel
!           Check if elements are completely outside searchbox
            inner: do kk = 1,nsd
!              Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).gt.sb%box(ii)%dim(2*kk  ))
     2            cycle outer
               if (elbox(2*kk  ,jj).lt.sb%box(ii)%dim(2*kk-1))
     2            cycle outer
            enddo inner

            sbel(cnt2) = jj
            cnt2=cnt2+1
         enddo outer
         ALLOCATE(sb%box(ii)%els(cnt2-1))
         sb%box(ii)%els=sbel(1:cnt2-1)

!        Now check face elements
!        Loop over faces
         outer2: do jj = 1,msh%nFa
         sbelf = 0
         cnt2  = 1
            middle: do kk = 1,msh%fa(jj)%nEl
                  inner2: do ll = 1,nsd
!                 Check if face elements are completely outside searchboxes
                  IF (MAXVAL(facels(jj)%elbox(2*ll  ,:)).lt.
     2            sb%box(ii)%dim(2*ll-1)) cycle outer2

                  IF (MINVAL(facels(jj)%elbox(2*ll-1,:)).gt.
     2            sb%box(ii)%dim(2*ll  )) cycle outer2

!                 Cycle if min value facels .gt. max value searchbox & vice-verse
                  if (facels(jj)%elbox(2*ll-1,kk).gt.
     2                  sb%box(ii)%dim(2*ll  )) cycle middle
                  if (facels(jj)%elbox(2*ll  ,kk).lt.
     2                  sb%box(ii)%dim(2*ll-1)) cycle middle
                  enddo inner2

                  sbelf(cnt2) = kk
                  cnt2 = cnt2 + 1
            enddo middle
            ALLOCATE(sb%box(ii)%fa(jj)%els(cnt2-1))
            sb%box(ii)%fa(jj)%els = sbelf(1:cnt2-1)
         enddo outer2
      end do
      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSb
!--------------------------------------------------------------------
!     Returns the ID of the searchboxes that contains point x
      PURE FUNCTION idSB(sb,x) RESULT(iSb)
      IMPLICIT NONE
      CLASS(sbType), INTENT(IN) :: sb
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER iSb(2**nsd)
      REAL(KIND=8) :: xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Set domain back to zero
      xzero(1) = x(1) - minval(sb%box(:)%dim(1))
      xzero(2) = x(2) - minval(sb%box(:)%dim(3))
      xzero(3) = x(3) - minval(sb%box(:)%dim(5))

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
      END FUNCTION idSB
!--------------------------------------------------------------------
!     Finds the element ID the particle is in. Also returns shape functions
      FUNCTION shapeFPrt(prt, ip, msh) RESULT(N)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN),TARGET :: prt
      INTEGER, INTENT(IN) :: ip
      TYPE(mshType), INTENT(IN) :: msh
      REAL(KIND=8) :: N(msh%eNoN)
      TYPE(pRawType), POINTER :: p
      TYPE(boxType),  POINTER :: b
      INTEGER :: ii, ind

      p => prt%dat(ip)
      b => prt%sb%box(MINVAL(p%sbID, MASK = p%sbID.gt.0))

      do ii=1,size(b%els)+1
            
            IF (ii.eq.1) THEN
                  ind = p%eIDo
            ELSE
                  ind = b%els(ii-1)
            END IF
            N = msh%nAtx(p%x,msh%x(:,msh%IEN(:,ind)))

            ! Checking if all shape functions are positive
            IF (ALL(N.ge.-1.0D-7)) then
                  p%eID=ind
                  p%N = N
                  EXIT
            END IF
         
      end do

      ! If it loops through everything and doesn't yield a positive shape function,
      ! the particle is outside the domain

      RETURN
      END FUNCTION shapeFPrt
!--------------------------------------------------------------------
      SUBROUTINE seedPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      
      INTEGER ip

      TYPE(pRawType) p

      CALL RSEED(cm%id())
      DO ip=1, prt%n
           p%u    = 0D0
           p%pID  = INT(ip,8)
           CALL RANDOM_NUMBER(p%x(:))
           IF (nsd .EQ. 2) p%x(3) = 0D0
           !p%x = (p%x - (/0.5D0,0.5D0,0D0/))*2D0
           p%x(1) = 0D0
           p%x(2) = 0.001D0
           p%x(3) = 0.31D0/ip
           prt%dat(ip) = p
      END DO


      RETURN
      END SUBROUTINE seedPrt

!--------------------------------------------------------------------
!     Gets the acceleration due to drag on a single particle
      FUNCTION dragPrt(prt, ip) RESULT(apd)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER,INTENT(IN) :: ip
      REAL(KIND=8) :: taupo      
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
      do ii=1,nsd
         do jj=1,msh%eNoN
            fvel(ii) = fvel(ii) + u%v(ii,msh%IEN(jj,p%eID))*p%N(jj)
         end do
      end do
      fvel = 0

      ! Relative velocity
      relvel = fvel - p%u
      ! Relative velocity magnitude
      magud = SUM(relvel**2D0)**0.5D0
      ! Reynolds Number
      Rep = dp*magud*rhoF/mu
      ! Schiller-Neumann (finite Re) correction
      fSN = 1D0 !+ 0.15D0*Rep**0.687D0   !!!
      ! Stokes corrected drag force
      apD = fSN/taup*relvel

      END FUNCTION dragPrt
!--------------------------------------------------------------------
!     Detects and collisions and returns pairs and times
      SUBROUTINE findcollPrt(prt,id1,id2,lM)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      CLASS(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: id1, id2
      TYPE(pRawType), POINTER :: p1,p2

!     Calculating distance coefficient
      REAL(KIND=8) :: a, b, c, d, e, f, qa, qb,qc,zeros(2),tcr,dp
      REAL(KIND=8) :: Np1(nsd+1), Np2(nsd+1)

      p1 => prt%dat(id1)
      p2 => prt%dat(id2)
      dp  = prt%mat%D()

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
      Np1 = prt%shapeF(id1, lM)
      Np2 = prt%shapeF(id2, lM)

      p1%x = p1%xo
      p1%x = p1%xo

      IF (ANY(Np1.lt.-1D-7) .or. ANY(Np2.lt.-1D-7)) THEN
            RETURN
      END IF

      p1%ti = tcr
      p2%ti = tcr

      p1%collided = .true.
      p2%collided = .true.

      prt%collcnt = prt%collcnt+1

      END SUBROUTINE findcollPrt

!--------------------------------------------------------------------
      SUBROUTINE collidePrt(prt,id1,id2)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: id1, id2
      TYPE(pRawType), POINTER :: p1,p2

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

      p1%collided = .false.
      p2%collided = .false.


      END SUBROUTINE collidePrt
!--------------------------------------------------------------------
      SUBROUTINE findwlPrt(prt,idp,msh)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      CLASS(mshType), INTENT(IN) :: msh
      TYPE(pRawType), POINTER :: p
      TYPE(boxType),  POINTER :: b
      REAL(KIND=8) :: Jac, xXi(nsd,nsd), Am(nsd,nsd), x1(nsd), tc
      REAL(KIND=8) :: N(msh%eNoN),xi(nsd),Bm(nsd), xc(nsd) 
      INTEGER :: ii, jj, a,gEl, faceNS(1,msh%fa(1)%eNoN), cnt,kk,ll
     2 , faceN(msh%fa(1)%eNoN),facev
      REAL(KIND=8) s, t, mx, my, ux, uy, uz, lx, ly, lz, iJac,
     2 xl(nsd,msh%eNoN)

      p => prt%dat(idp)
      b => prt%sb%box(MINVAL(p%sbID, MASK = p%sbID.gt.0))

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
            p%x = (/0D0,0D0,29.5D0/)!msh%x(:,msh%fa(ii)%IEN(1,rndi))

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
      TYPE(matType), POINTER :: mat
      TYPE(mshType), POINTER :: msh
      TYPE(pRawType), POINTER :: p, tp
      TYPE(prtType), TARGET :: tmpprt
      TYPE(sbType), POINTER :: sb
      TYPE(VarType),POINTER :: u
      REAL(KIND=8) rhoF,
     2   mu, mp, g(nsd), dp,
     3   rhoP, N(nsd+1), Ntmp(nsd+1), tt(4)
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd)
      REAL(KIND=8) :: prtxpred(nsd), pvelpred(nsd), tmpwr(nsd),mom
      INTEGER :: a, Ac

!     Gravity
!!!!! Find where grav actually is? (maybe mat%body forces)
      g=0D0
      g(3)=-1D0
      if (prt%ptr(idp) .eq. 2) g=-g
      tmpprt = prt
      msh => prt%dmn%msh(1)
      p  => prt%dat(idp)
      sb => prt%sb
      tp => tmpprt%dat(1)
      u => prt%Uns
      tmpprt%sb = sb
      rhoF  = prt%mns%rho()
      rhoP  = prt%mat%rho()
      mu    = prt%mns%mu()
      dp    = prt%mat%D()
      mp = pi*rhoP/6D0*dp**3D0

      p%eIDo = p%eID
      p%sbIDo= p%sbID


      IF (p%wall) THEN
!           Find which searchbox particle is in
            p%sbID = sb%id(p%x)
!           Get shape functions/element of particle
            N = prt%shapeF(idp, msh)
      END IF

!     Get drag acceleration
      apd = prt%drag(idp)
!     Total acceleration (just drag and buoyancy now)
      apT = apd + g*(1D0 - rhoF/rhoP)
!     2nd order advance (Heun's Method)
!     Predictor
      pvelpred = p%u + p%remdt*apT
      prtxpred = p%x + p%remdt*p%u
      tp%u  = pvelpred
      tp%x  = prtxpred

!     Find which searchbox prediction is in
      tp%sbID = sb%id(tp%x)   !!!!!!!!!!!!! For some reason this changes apT???? Does not make any sense to me at all

!     Get shape functions/element of prediction
      Ntmp = tmpprt%shapeF(1, msh)
!     Check if predictor OOB
      IF (ANY(Ntmp.le.-1D-7)) THEN
!     This should advance to the edge of the wall, and then change the velocitry, as well as giving remdt
            CALL prt%findwl(idp,msh)
            CALL prt%wall(idp,msh)
            p%x = p%x + p%remdt*p%u
            p%wall = .TRUE.
            RETURN
      END IF
!     Total acceleration (just drag and buoyancy now) ! again, because above changes this for some reason
      apT = apd + g*(1D0 - rhoF/rhoP)

!     Get drag acceleration of predicted particle
      apdpred = tmpprt%drag(1)
      apTpred = apdpred + g*(1D0 - rhoF/rhoP)
!     Corrector
      p%u = p%u + 0.5D0*p%remdt*(apT+apTpred)
      p%x = p%remdt*p%u + p%x

!     Send drag to fluid
      DO a=1,msh%eNoN
            Ac = msh%IEN(a,p%eID)
!            prt%twc%v(:,Ac) = prt%twc%v(:,Ac) +
!     2      0.5D0*(apd + apdpred)*mP/rhoF/prt%wV(Ac)*p%N(a)
      END DO

      IF (time.lt.0.04) prt%twc%v(1,1) = 1D0

      IF (prt%itr .EQ. 0) THEN
            tmpwr = 0.5*(apd+apdpred)
            write(88,*) sqrt(tmpwr(1)**2+tmpwr(2)**2+tmpwr(3)**2)*mp
            !print *, sqrt(tmpwr(1)**2+tmpwr(2)**2+tmpwr(3)**2)*mp
            !mom =prt%dmn%msh(1)%integ(u%v, 3)
            !print *, mom
            !call sleep(1)
      END IF

!     Check if particle went out of bounds
      p%sbID = sb%id(p%x)
      N = prt%shapeF(idp, msh)

      IF (ANY(N .le. -1D-7)) THEN
            p%x = p%xo
            CALL prt%findwl(idp,msh)
            CALL prt%wall(idp,msh)
            p%x = p%x + p%remdt*p%u
            p%wall = .TRUE.
            RETURN
      END IF

      p%wall = .FALSE.
      RETURN

      END SUBROUTINE advPrt
!--------------------------------------------------------------------
      SUBROUTINE solvePrt(eq)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT):: eq
      INTEGER ip, a, e, eNoN, Ac,i ,j, k,l, subit
      TYPE(mshType), POINTER :: lM
!     Particle/fluid Parameters
      REAL(KIND=8):: dtp,maxdtp,sbdt(nsd),dp,taup,rhoP,mu,tim,
     2 P1, P2

      lM => eq%dmn%msh(1)
      rhoP  = eq%mat%rho()
      mu    = eq%mns%mu()
      dp    = eq%mat%D()

!     Reset twc force to zero
      eq%twc%v(:,:) = 0D0

!     Reset collision counter
      eq%collcnt = 0

!     Particle relaxation time
      taup = rhoP*dp**2D0/mu/18D0

!!!!! get time step broken down to be dictated by either fastest velocity(in terms of eq%sb's traveled), relax time, overall solver dt
!!! idea: do one more loop through all partsicles above this one and get time step for each one, then take minimum
!!! Right now: I'm just going to take min of relaxation time and overall, but will need to make it eq%sb so I can only check neighboring eq%sb's

!     Appropriate time step
      dtp = MIN(taup,dt)
!     Subiterations
      subit = FLOOR(dt/dtp)
      dtp = dt/subit
      eq%dt = dtp

      DO l = 1,subit

      DO i=1,eq%n

!     If it's the first iteration, update last velocity, position
            IF (eq%itr .EQ. 0) THEN
                  eq%dat(i)%xo = eq%dat(i)%x
                  eq%dat(i)%uo = eq%dat(i)%u
            ENDIF

!     Set position and velocity to old variables in preparation for iteration with INS
            eq%dat(i)%x = eq%dat(i)%xo
            eq%dat(i)%u = eq%dat(i)%uo

!     Set initial advancing time step to solver time step
            eq%dat(i)%remdt = dtp
      ENDDO

      DO i = 1,eq%n
      !     Collisions
            DO j=1,eq%n
!     Check if the particle collides with any other particles and hasn't collided. Advance to collision location
            IF ((i.ne.j).and.(.not.(eq%dat(i)%collided))) THEN
                  CALL eq%findcoll(i,j,lM)
                  IF (eq%dat(i)%collided) CALL eq%collide(i,j)
            END IF
            ENDDO
      ENDDO

      DO i = 1,eq%n
            CALL eq%adv(i)
            IF (eq%itr .EQ. 0)  print *, eq%dat(i)%x(3),
     2            eq%dat(i)%u(3) ,taup    
!           Reset if particles have collided
            eq%dat(i)%collided = .false.
      END DO
            
      P1 = eq%dmn%msh(1)%integ(1,eq%Pns%s)
      P2 = eq%dmn%msh(1)%integ(2,eq%Pns%s)

      IF (eq%itr .EQ. 0) write(88,*) eq%dat(1)%u(3), !eq%dat(2)%u
     2 eq%dmn%msh(1)%integ(eq%Uns%v, 3)*1.2D0,
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


      !!!! Urgent fixes after wall:
      !!!! Add in so it only checks in same searchbox
      !!!! Mult collisions per step
      !!! Right now doesn't consider collisions after other collisions (wall or particle)