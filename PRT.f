      MODULE PRTMOD
      USE INSMOD
!     FSI is using PRTMOD
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: OLD=1, CUR=2, NEW=3, RHS=4, IMP=4

!     Size of the container is extended by this factor if it's
!     completely filled
      INTEGER, PARAMETER :: prtExtFac = 2

!     Face element boxes
      TYPE facelsboxType
            REAL(KIND=8), ALLOCATABLE :: elbox(:,:)
      END TYPE

!     Faces and their elements
      TYPE facelsType
      ! Elements of face
            INTEGER, ALLOCATABLE :: els(:)
      ! Type of face (1 = wall, 2 = inlet, 3 = outlet)
            INTEGER :: type = 0

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
         REAL(KIND=8) N(4)
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
      TYPE prtType
!     Created?
         LOGICAL :: crtd = .FALSE.
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
!     Material properties
         TYPE(matType), POINTER :: mat

      CONTAINS
!     Creates a new prt type
         PROCEDURE :: new => newPrt
!     Selects one particle
         PROCEDURE :: get => getFPrt
!     Removes a particle
         PROCEDURE :: rm => rmFPrt
!     Adds a particle to the structure
         PROCEDURE :: add => addFPrt
!     Increases the size of the container 
         PROCEDURE :: extend => extendPrt
!     Deletes the structure
         PROCEDURE :: free => freePrt
!     Returns shape function values at the location of a particle
         PROCEDURE :: shapeF => shapeFPrt
!     Seed the domai with particles
         PROCEDURE :: seed => seedPrt
!     Finds accelereation on a particle from drag
         PROCEDURE :: drag => dragPrt
!     Advance all particles one time step
         PROCEDURE :: solve => solvePrt
!     Detects and enacts collisions
         PROCEDURE :: collide => collidePrt
!     Finds which wall element the partice collides with
         PROCEDURE :: findwl => findwlPrt
!     Enacts collisions with walls
         PROCEDURE :: wall => wallPrt
!     Advances particle given time step
         PROCEDURE :: adv => advPrt
      END TYPE prtType

      CONTAINS

!####################################################################

      SUBROUTINE newPrt(prt, n)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      INTEGER, INTENT(IN), OPTIONAL :: n

      INTEGER i

      IF (prt%crtd) RETURN
      open(88,file='pos.txt')
      IF (PRESENT(n)) prt%maxN = n
      ALLOCATE(prt%dat(prt%maxN), prt%ptr(prt%maxN))
      DO i=1, prt%maxN
         prt%ptr(i) = i
      END DO
      prt%n = n
      prt%mat  => FIND_MAT('Particle')
      CALL prt%seed()
      prt%crtd = .TRUE.

      RETURN
      END SUBROUTINE newPrt
!--------------------------------------------------------------------
      FUNCTION getFPrt(d,i)
      IMPLICIT NONE
      CLASS(prtType), TARGET, INTENT(IN) :: d
      INTEGER, INTENT(IN) :: i
      TYPE(pRawType), POINTER :: getFPrt

      INTEGER ind

      IF (i.GT.d%n .OR. i.LT.1) THEN
         getFPrt => NULL()
         RETURN
      END IF

      ind = d%ptr(i)
      getFPrt => d%dat(ind)

      RETURN
      END FUNCTION getFPrt
!--------------------------------------------------------------------
      SUBROUTINE rmFPrt(d,i)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: d
      INTEGER, INTENT(IN) :: i

      INTEGER ind

      IF (.NOT.d%crtd) STOP "prt not created yet!"
      IF (i.GT.d%n .OR. i.LT.1) STOP "Out of range index in rmFPrt"

      ind        = d%ptr(i)
      d%ptr(i)   = d%ptr(d%n)
      d%ptr(d%n) = ind
      d%n        = d%n - 1

      RETURN
      END SUBROUTINE rmFPrt
!--------------------------------------------------------------------
      SUBROUTINE addFPrt(d,iDat,ind)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: d
      TYPE(pRawType), INTENT(IN) :: iDat

      INTEGER, INTENT(IN) :: ind

      IF (.NOT.d%crtd) STOP "prt not created yet!"
      !IF (d%n .EQ. d%maxN) CALL d%extend()

      !d%n        = d%n + 1
      !ind        = d%ptr(d%n)
      d%dat(ind) = iDat

      RETURN
      END SUBROUTINE addFPrt
!--------------------------------------------------------------------
      SUBROUTINE extendPrt(d, maxN)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: d
      INTEGER, INTENT(IN), OPTIONAL :: maxN

      INTEGER oldMax, newMax
      TYPE(prtType) dTmp

      IF (.NOT.d%crtd) STOP "prt not created yet!"
      
      oldMax = d%maxN
      newMax = prtExtFac*oldMax
      IF (PRESENT(maxN)) newMax = maxN
      IF (newMax .EQ. 0) newMax = prtExtFac

!     At least size must be increased by one
      IF (newMax .LE. oldMax) RETURN

      CALL dTmp%new(newMax)
      dTmp%n             = d%n
      dTmp%dat(1:oldMax) = d%dat
      dTmp%ptr(1:oldMax) = d%ptr

      CALL d%free()
      CALL d%new(newMax)
      d%n             = dTmp%n
      d%dat(1:oldMax) = dTmp%dat(1:oldMax)
      d%ptr(1:oldMax) = dTmp%ptr(1:oldMax)
      CALL dTmp%free()

      RETURN
      END SUBROUTINE extendPrt
!--------------------------------------------------------------------
      PURE SUBROUTINE freePrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt

      INTEGER iSb

      IF (.NOT.prt%crtd) RETURN
      
      prt%crtd = .FALSE.
      prt%n    = 0
      DEALLOCATE(prt%dat, prt%ptr)
      

      RETURN
      END SUBROUTINE freePrt
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

      split=(/2,2,2/)
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
      ! Tolerance
      eps = 0.5*diff
      ! Size of sb
      sb%step=diff/((sb%n+1)/2) + eps

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

!     Make boxes around face elements
!!!! This part seems right
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
      INTEGER :: xsteps(nsd), i(nsd)

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
      REAL(KIND=8) :: N(4)
      TYPE(mshType), INTENT(IN) :: msh

      INTEGER, PARAMETER :: TIME_POINT=CUR, iM=1

      LOGICAL NOX
      TYPE(pRawType), POINTER :: p
      TYPE(boxType),  POINTER :: b
      INTEGER :: ii,cnt,a, ind
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,nsd+1),
     2 Nxi(nsd,nsd+1),prntx(nsd)
      cnt=1

      Nxi(1,1) =  1D0
      Nxi(2,1) =  0D0
      Nxi(3,1) =  0D0
      Nxi(1,2) =  0D0
      Nxi(2,2) =  1D0
      Nxi(3,2) =  0D0
      Nxi(1,3) =  0D0
      Nxi(2,3) =  0D0
      Nxi(3,3) =  1D0
      Nxi(1,4) = -1D0
      Nxi(2,4) = -1D0
      Nxi(3,4) = -1D0

      IF (msh%eType.NE.eType_TET) 
     2   io%e = "shapeFPrt only defined for tet elements"

      p => prt%dat(ip)
      b => prt%sb%box(p%sbID(1))

      do ii=1,size(b%els)+1

            ind = b%els(ii-1)
            if (ii.eq.1) ind = p%eIDo
            xXi = 0D0

      IF (nsd .EQ. 2) THEN
      !
      ! 2D not done
      !
!         DO a=1, eNoN
!            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
!            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
!         END DO
!
!         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)
!
!         xiX(1,1) =  xXi(2,2)/Jac
!         xiX(1,2) = -xXi(1,2)/Jac
!         xiX(2,1) = -xXi(2,1)/Jac
!         xiX(2,2) =  xXi(1,1)/Jac
!
!         
!         DO a=1, eNoN
!            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
!            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
!         END DO
!
      ELSE
      ! 3D case

      ! Setting up matrix to be inverted
         DO a=1, msh%eNoN
            xXi(:,1) = xXi(:,1) +
     2        msh%x(:,msh%IEN(a,ind))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) +
     2        msh%x(:,msh%IEN(a,ind))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + 
     2        msh%x(:,msh%IEN(a,ind))*Nxi(3,a)
         END DO
         
      ! Inverting matrix
         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2       + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3       + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4       - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5       - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6       - xXi(1,3)*xXi(2,2)*xXi(3,1)

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac
         
      ! Finding particle coordinates in parent coordinate system
         DO a=1, nsd
            prntx(a) = xiX(a,1)*(p%x(1)- msh%x(1,msh%IEN(4,ind)))
     2               + xiX(a,2)*(p%x(2)- msh%x(2,msh%IEN(4,ind)))
     3               + xiX(a,3)*(p%x(3)- msh%x(3,msh%IEN(4,ind)))
         END DO
      END IF

      ! Finding shape function values at particle coordinates
      N(1) = prntx(1)
      N(2) = prntx(2)
      N(3) = prntx(3)
      N(4) = 1 - prntx(1) - prntx(2) - prntx(3)

      ! Checking if all shape functions are positive
      IF (ALL(N.ge.0D0)) then
         p%eID=ind
         p%N = N
         EXIT
      END IF
         
      end do

      ! If it loops through everything and doesn't yield a positive shape function,
      ! the particle is outside the domain
      !if (p%eID.eq.0) io%e = 'outside domain'

      RETURN
      END FUNCTION shapeFPrt
!--------------------------------------------------------------------
      SUBROUTINE seedPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      
      INTEGER ip, iDis

      TYPE(pRawType) p

      CALL RSEED(cm%id())
      DO ip=1, prt%n
           p%u    = 0D0
           p%pID  = INT(ip,8)
           CALL RANDOM_NUMBER(p%x(:))
           IF (nsd .EQ. 2) p%x(3) = 0D0
           !p%x = (p%x - (/0.5D0,0.5D0,0D0/))*2D0
           !CALL prt%add(p)
           p%x(1) = 0D0
           p%x(2) = 0D0
           p%x(3) = 5D0/ip
           prt%dat(ip) = p
      END DO


      RETURN
      END SUBROUTINE seedPrt
!--------------------------------------------------------------------
      SUBROUTINE PRTTRANS
      IMPLICIT NONE 

      INTEGER e, ip, a, Ac, eNoN
      REAL(KIND=8) N(4), up(3)
      TYPE(prtType) :: prt
      TYPE(pRawType), POINTER :: p
      TYPE(mshType) :: msh
      ! wrong here vv
      REAL(KIND=8) :: YN(3,1)
 
      eNoN = msh%eNoN
!     Initializing particles
      CALL prt%new()
      CALL prt%seed()
      CALL prt%sb%new(msh)



      DO ip=1, prt%n
         p => prt%get(ip)
         !e = prt%shapeF(ip,msh)
         up = 0D0
         DO a=1, eNoN
            Ac = msh%IEN(a,e)
            ! Don't know what Yn is
            up = up + Yn(1:3,Ac)*N(a)
         END DO
      END DO

      call prt%free()

      RETURN
      END SUBROUTINE PRTTRANS
!--------------------------------------------------------------------
!     Gets the acceleration due to drag on a single particle
      FUNCTION dragPrt(prt, ip, ns) RESULT(apd)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER,INTENT(IN) :: ip
      CLASS(eqType), INTENT(IN) :: ns
      REAL(KIND=8) :: taupo      
      TYPE(gVarType),POINTER :: u
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
      mat => ns%mat
      msh => ns%dmn%msh(1)
      u => ns%var(1)

!     Fluid parameters
      rhoF = mat%rho()
      mu  = mat%mu()

!     Particle relaxation time
      taup = rhoP*dp**2D0/mu/18D0

!     Interpolate velocity at particle point
      fvel=0D0
      do ii=1,nsd
         do jj=1,msh%eNoN
            fvel(ii) = fvel(ii) + u%v(ii,msh%IEN(jj,p%eID))*p%N(jj)
         end do
      end do
      !fvel = 0
      !fvel(3) = -3D0!!!!!!!!!!!!!!
      !fvel(1) = fvel(1)+1
      !if (ip.eq.2) fvel(3) = 3D0!!!!!!!!

      ! Relative velocity
      relvel = fvel - p%u
      ! Relative velocity magnitude
      magud = SUM(relvel**2D0)**0.5D0
      ! Reynolds Number
      Rep = dp*magud*rhoP/mu
      ! Schiller-Neumann (finite Re) correction
      fSN = 1D0 + 0.15D0*Rep**0.687D0
      ! Stokes corrected drag force
      apD = fSN/taup*relvel

      END FUNCTION dragPrt
!--------------------------------------------------------------------
!     Detects and enacts collisions
      SUBROUTINE collidePrt(prt,id1,id2,dtp,lM)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      CLASS(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: id1, id2
      REAL(KIND=8), INTENT(IN) :: dtp
      TYPE(pRawType), POINTER :: p1,p2

!     Calculating distance coefficient
      REAL(KIND=8) :: a, b, c, d, e, f, qa, qb, qc, zeros(2), tcr
      REAL(KIND=8) :: n1(nsd), n2(nsd), t1(nsd), t2(nsd)
      REAL(KIND=8) :: vpar1, vpar2, vperp1, vperp2, dp, mp, rho, k
!     Coefficients to make calculating parallel/perp vel easier
      REAL(KIND=8) :: pa, pb, Np1(nsd+1), Np2(nsd+1) , temp(nsd)

      p1 => prt%dat(id1)
      p2 => prt%dat(id2)
      dp  = prt%mat%D()
      rho = prt%mat%rho()
      k   = prt%mat%krest()
      mp = pi*rho/6D0*dp**3D0

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

      ! Exit function if collision won't occur during timestep
      if (tcr.gt.dtp) RETURN

      ! particle locations at point of collision
      p1%xc = p1%u*tcr + p1%x
      p2%xc = p2%u*tcr + p2%x
      p1%x  = p1%xc
      p2%x  = p2%xc

!     Check if the particle is outside the domain
      Np1 = prt%shapeF(id1, lM)
      Np2 = prt%shapeF(id2, lM)
      IF (ANY(Np1.lt.0) .or. ANY(Np2.lt.0)) THEN
            p1%x = p1%xo
            p2%x = p2%xo
            RETURN
      END IF

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

      p1%u = vpar1*n1 + vperp1*t1
      p2%u = vpar2*n2 + vperp2*t2
            
      ! particle velocities at point of collision
      p1%uc = p1%u
      p2%uc = p2%u

      p1%collided = .true.
      p2%collided = .true.

      p1%remdt = p1%remdt-tcr
      p2%remdt = p2%remdt-tcr

      END SUBROUTINE collidePrt
!--------------------------------------------------------------------
      SUBROUTINE findwlPrt(prt,idp,msh)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      CLASS(mshType), INTENT(IN) :: msh
      TYPE(pRawType), POINTER :: p
      TYPE(boxType),  POINTER :: b
      REAL(KIND=8) :: Jac, xXi(nsd,nsd), xiX(nsd,nsd), Nxi(nsd,nsd)
      REAL(KIND=8) :: prntx(nsd), N(nsd)
      INTEGER :: ii, jj, a

      p => prt%dat(idp)
      b => prt%sb%box(p%sbID(1))

      Nxi(1,1) =  1D0
      Nxi(2,1) =  0D0
      Nxi(1,2) =  0D0
      Nxi(2,2) =  1D0
      Nxi(1,3) = -1D0
      Nxi(2,3) = -1D0
      p%faID = 0

      faceloop: DO ii=1,msh%nFa
      DO jj=1,size(b%fa(ii)%els)

            xXi = 0D0
      !     Setting up matrix for inversion
            DO a=1, msh%eNoN-1
                  xXi(:,1) = xXi(:,1) +
     2        msh%x(:,msh%fa(ii)%IEN(a,b%fa(ii)%els(jj)))*Nxi(1,a)
                  xXi(:,2) = xXi(:,2) +
     2        msh%x(:,msh%fa(ii)%IEN(a,b%fa(ii)%els(jj)))*Nxi(2,a)
            ENDDO
            xXi(:,3) = -p%u

                  ! Inverting matrix
            Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2       + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3       + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4       - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5       - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6       - xXi(1,3)*xXi(2,2)*xXi(3,1)
    
            xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
            xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
            xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
            xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
            xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
            xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
            xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
            xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
            xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac

            DO a=1, nsd
            prntx(a) = xiX(a,1)*(p%x(1) -
     2       msh%x(1,msh%fa(ii)%IEN(3,b%fa(ii)%els(jj))))
     3               + xiX(a,2)*(p%x(2) -
     4       msh%x(2,msh%fa(ii)%IEN(3,b%fa(ii)%els(jj))))
     5               + xiX(a,3)*(p%x(3) -
     6       msh%x(3,msh%fa(ii)%IEN(3,b%fa(ii)%els(jj))))
            END DO

            N(1) = prntx(1)
            N(2) = prntx(2)
            N(3) = 1 - prntx(1) - prntx(2)

            IF (ALL(N.ge.0D0).and. prntx(3).gt.0
     2      .and. prntx(3).lt.prt%dt) THEN
                  p%faID(1) = ii
                  p%faID(2) = b%fa(ii)%els(jj)
                  p%ti = prntx(3)
                  EXIT faceloop
            ENDIF

      ENDDO
      ENDDO faceloop
      !io%e = "Wrong searchbox"
      
      END SUBROUTINE findwlPrt
!--------------------------------------------------------------------
      SUBROUTINE wallPrt(prt, idp, msh)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      CLASS(mshType), INTENT(IN) :: msh
      TYPE(pRawType), POINTER :: p
      INTEGER :: ii, rndi
      REAL(KIND=8) :: dp, rho, k, nV(nsd), tV(nsd),
     2 a(nsd), b(nsd), vpar, vperp, temp(nsd), rnd

      p   =>prt%dat(idp)
      dp  = prt%mat%D()
      rho = prt%mat%rho()
      k   = prt%mat%krest()

!     Advance to collision location
      p%xc = p%u*p%ti + p%x
      p%x  = p%xc
      p%remdt = p%remdt - p%ti

      IF (msh%fa(p%faID(1))%typ .EQ. 3) THEN
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
!           If inlet...
            IF (msh%fa(ii)%typ .EQ. 1) THEN
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
      SUBROUTINE advPrt(prt, idp, ns)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: idp
      CLASS(eqType), INTENT(IN) :: ns
      !TYPE(insType), POINTER :: ins
      TYPE(matType), POINTER :: mat
      TYPE(mshType), POINTER :: msh
      TYPE(pRawType), POINTER :: p, tp
      TYPE(prtType), TARGET :: tmpprt
      TYPE(sbType), POINTER :: sb
      TYPE(gVarType),POINTER :: u
      REAL(KIND=8) rhoF,
     2   mu, mp, g(nsd), dp,
     3   rhoP, N(nsd+1), Ntmp(nsd+1)
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd)
      REAL(KIND=8) :: prtxpred(nsd), pvelpred(nsd)
      INTEGER :: a

!     Gravity
!!!!! Find where grav actually is? (maybe mat%body forces)
      g(3)=10D0
      !ins => ns%s

      tmpprt = prt
      mat => ns%mat
      msh => ns%dmn%msh(1)
      p  => prt%dat(idp)
      sb => prt%sb
      tp => tmpprt%dat(1)
      u => ns%var(1)
      tmpprt%sb = sb
      rhoF  = ns%mat%rho()
      rhoP  = prt%mat%rho()
      mu    = ns%mat%mu()
      dp    = prt%mat%D()
      mp = pi*rhoP/6D0*dp**3D0

!     Set last coordinates
      p%xo = p%x
      p%uo = p%u
      p%eIDo = p%eID
      p%sbIDo= p%sbID

 1    CONTINUE

      IF (p%wall) THEN
!           Find which searchbox particle is in
            p%sbID = sb%id(p%x)
!           Get shape functions/element of particle
            N = prt%shapeF(idp, msh)
      END IF

!     Get drag acceleration
      apd = prt%drag(idp,ns)
!     Total acceleration (just drag and buoyancy now)
      apT = apd + g*(1D0 - rhoF/rhoP)
!     2nd order advance (Heun's Method)
!     Predictor
      pvelpred = p%u + p%remdt*apT
      prtxpred = p%x + p%remdt*p%u
      tp%u  = pvelpred
      !if (idp.eq.2) tp%u = -pvelpred!!!!!!!!!!!!!!
      tp%x  = prtxpred

!     Find which searchbox prediction is in
      tp%sbID = sb%id(tp%x)

!     Get shape functions/element of prediction
      Ntmp = tmpprt%shapeF(1, msh)

!     Check if predictor OOB
      IF (ANY(Ntmp.lt.0)) THEN
!     This should advance to the edge of the wall, and then change the velocitry, as well as giving remdt
            CALL prt%findwl(idp,msh)
            CALL prt%wall(idp,msh)
            p%x = p%x + p%remdt*p%u
            p%wall = .TRUE.
            RETURN
      END IF

!     Get drag acceleration of predicted particle
      apdpred = tmpprt%drag(1,ns)
      !if (idp.eq.2) apdpred = -apdpred!!!!!!!!!!!!!!!
      apTpred = apdpred + g*(1D0 - rhoF/rhoP)
!     Corrector
      p%u = p%u + 0.5D0*p%remdt*(apT+apTpred)
      p%x = p%remdt*p%u + p%x

      DO a=1,msh%eNoN
            u%OC%v(:,msh%IEN(a,p%eID)) = 
     3      -0.5*(apd + apdpred)*rhoP/rhoF*p%N(a)
      END DO
      !print *, ns%var(1)%OC%v(1,p%eID)
      !print *, u%v(:,msh%IEN(3,p%eIDo))

!     Check if particle went out of bounds
      p%sbID = sb%id(p%x)
      N = prt%shapeF(idp, msh)

      IF (ANY(N .lt. 0)) THEN
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
      SUBROUTINE solvePrt(prt, ns)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
      CLASS(eqType), INTENT(IN) :: ns
      INTEGER ip, a, e, eNoN, Ac,i ,j, k,l, subit
      TYPE(sbType), POINTER :: sb
      TYPE(mshType), POINTER :: lM
      TYPE(pRawType), POINTER :: p(:)
      TYPE(matType), POINTER :: mat
!     Particle/fluid Parameters
      REAL(KIND=8) :: dtp,maxdtp,sbdt(nsd), dp, taup,rhoP,mu

!     Initialize if haven't yet
      IF(.NOT.prt%crtd) THEN
            CALL prt%new(1)
      END IF      

      lM => ns%dmn%msh(1)
      p  => prt%dat
      sb => prt%sb
      mat => ns%mat
      rhoP  = prt%mat%rho()
      mu    = ns%mat%mu()
      dp    = prt%mat%D()

!     Particle relaxation time
      taup = rhoP*dp**2D0/mu/18D0

      IF (.NOT.(sb%crtd)) THEN
         CALL sb%new(lM)
      END IF

!!!!! get time step broken down to be dictated by either fastest velocity(in terms of sb's traveled), relax time, overall solver dt
!!! idea: do one more loop through all particles above this one and get time step for each one, then take minimum
!!! Right now: I'm just going to take min of relaxation time and overall, but will need to make it sb so I can only check neighboring sb's

!     Appropriate time step
      dtp = MIN(taup,dt)
!     Subiterations
      subit = FLOOR(dt/dtp)
      dtp = dt/subit
      prt%dt = dtp

      DO l = 1,subit

      !!! Do all collisions first, and only advance to point of collision
      !!! 2nd order: get vel at t impact, use that, but
      !!! Do same for wall collisions if it is out before coll
      DO i=1,prt%n
      p(i)%remdt = dtp
!        Collisions
         DO j=1,prt%n
!        Check if the particle collides with any other particles and hasn't collided. Advance if so.
            IF ((i.ne.j).and.(.not.(p(i)%collided))) THEN
                CALL prt%collide(i,j,dtp,lM)
            END IF
         ENDDO
      ENDDO

      DO i = 1,prt%n
            CALL prt%adv(i, ns)
            print *, p(i)%x
!        Reset if particles have collided
         p(i)%collided = .false.
      END DO

      write(88,*) p(1)%x!, p(2)%x

      ENDDO

      RETURN
      END SUBROUTINE solvePrt
      
      END MODULE PRTMOD


      !!!! Urgent fixes after wall:
      !!!! Add in so it only checks in same searchbox
      !!!! Velocity seems weird
      !!!! check if it is putting momentum back in

      !!! Right now doesn't consider collisions after other collisions (wall or particle)