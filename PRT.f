      MODULE PRTMOD
      USE EQMOD
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: OLD=1, CUR=2, NEW=3, RHS=4, IMP=4

!     Size of the container is extended by this factor if it's
!     completely filled
      INTEGER, PARAMETER :: prtExtFac = 2

!     Individual box
      TYPE boxType
!     Box dimensions (minx,maxx,miny,maxy,minz,maxz)
            REAL(KIND=8), ALLOCATABLE :: dim(:)
!     Elements contained in searchbox
            INTEGER, ALLOCATABLE :: els(:)
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
         INTEGER(KIND=8) :: eID=0
!     Particle ID
         INTEGER(KIND=8) :: pID=0
!     Searchbox ID
         INTEGER(KIND=8) :: sbID(8) = 0
!     Position
         REAL(KIND=8) x(3)
!     Velocity
         REAL(KIND=8) u(3)
!     Shape functions in current element
         REAL(KIND=8) N(4)
!     Has this particle collided during the current time step?
         LOGICAL :: collided = .FALSE.
!     Remaining time in timestep after collision
         REAL(KIND=8) :: remdt
      END TYPE pRawType

!     Collection of particles
      TYPE prtType
!     Created?
         LOGICAL :: crtd = .FALSE.
!     Maximum capacity of container
         INTEGER :: maxN = 1024
!     Current size
         INTEGER :: n = 0
!     Search boxes to find element hosting particles
         TYPE(sbType) :: sb
!     Array pointing to full and empty entries
         INTEGER, ALLOCATABLE :: ptr(:)
!     Data
         TYPE(pRawType), ALLOCATABLE :: dat(:)
!     Material properties
         TYPE(matType) :: mat

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
      END TYPE prtType

      CONTAINS

!####################################################################

      SUBROUTINE newPrt(prt, n)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: prt
      INTEGER, INTENT(IN), OPTIONAL :: n

      INTEGER i

      IF (prt%crtd) RETURN

      IF (PRESENT(n)) prt%maxN = n
      ALLOCATE(prt%dat(prt%maxN), prt%ptr(prt%maxN))
      DO i=1, prt%maxN
         prt%ptr(i) = i
      END DO
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
      SUBROUTINE addFPrt(d,iDat)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT) :: d
      TYPE(pRawType), INTENT(IN) :: iDat

      INTEGER ind

      IF (.NOT.d%crtd) STOP "prt not created yet!"
      IF (d%n .EQ. d%maxN) CALL d%extend()

      d%n        = d%n + 1
      ind        = d%ptr(d%n)
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
      INTEGER :: ii,jj,cnt2,kk
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), elbox(2*nsd,msh%nEl)
      INTEGER, ALLOCATABLE :: sbel(:)
      INTEGER :: split(nsd)
      IF (sb%crtd) RETURN

 !!     ! Will need to replace eventually with something that targets element size
      
      split=(/10,10,10/)
      sb%n=split

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%box(sb%n(1)*sb%n(2)*sb%n(3)))
      ALLOCATE(sbel(msh%nEl))


      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(sb%n(3)*sb%n(2)),seq2(sb%n(3)*sb%n(1))
     2   ,seq3(sb%n(2)*sb%n(1)))

      ! Domain ranges
      diff(1)=MAXVAL(msh%x(1,:))-MINVAL(msh%x(1,:))
      diff(2)=MAXVAL(msh%x(2,:))-MINVAL(msh%x(2,:))
      diff(3)=MAXVAL(msh%x(3,:))-MINVAL(msh%x(3,:))
      ! Size of sb
      sb%step=diff/((split+1)/2)

      seq1=(/(ii, ii=0, split(2)*split(3)-1, 1)/)*split(1)+1
      cnt2=0
      do ii=1,split(1)*split(3)
            seq2(ii)=ii+cnt2*(split(2)-1)*split(1)
            if (MOD(ii,split(1)).eq.0) cnt2=cnt2+1
      end do
      seq3=(/(ii, ii=0, split(1)*split(2)-1, 1)/)+1

      ! Allocating sb, such that they overlap by 50%
      ! Direction 1
      do ii=1,sb%n(1)
         sb%box(seq1+ii-1)%dim(1) = MINVAL(msh%x(1,:)) 
     2       + sb%step(1)*(ii-1)/2
      end do

      ! Direction 2
      do ii=1,sb%n(2)
         sb%box(seq2+(ii-1)*split(1))%dim(3) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)/2
      end do

      ! Direction 3
      do ii=1,sb%n(3)
         sb%box(seq3+(ii-1)*split(1)*split(2))%dim(5)=
     2   MINVAL(msh%x(3,:)) + sb%step(3)*(ii-1)/2
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

      do ii=1,sb%n(1)*sb%n(2)*sb%n(3)
         cnt2=1
         sbel=0
         outer: do jj=1,msh%Nel
            ! Check if elements are completely outside searchbox
            inner: do kk=1,nsd
               ! Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).gt.sb%box(ii)%dim(2*kk  )) cycle
               if (elbox(2*kk  ,jj).lt.sb%box(ii)%dim(2*kk-1)) cycle
            enddo inner

            sbel(cnt2) = jj
            cnt2=cnt2+1
         enddo outer
         ALLOCATE(sb%box(ii)%dim(2*nsd))
         ALLOCATE(sb%box(ii)%els(cnt2-1))
         sb%box(ii)%els=sbel(1:cnt2-1)
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
      xsteps = FLOOR(xzero/sb%step)
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
      FUNCTION shapeFPrt(prt, ip, msh) RESULT(e)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN),TARGET :: prt
      INTEGER, INTENT(IN) :: ip
      REAL(KIND=8) :: N(4)
      TYPE(mshType), INTENT(IN) :: msh
      INTEGER e

      INTEGER, PARAMETER :: TIME_POINT=CUR, iM=1

      LOGICAL NOX
      INTEGER  iSb(2**nsd)
      TYPE(pRawType), POINTER :: p
      TYPE(boxType),  POINTER :: b
      INTEGER :: ii,cnt,a
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,nsd+1),
     2 Nxi(nsd,nsd+1),prntx(nsd)
      cnt=1

      IF (msh%eType.NE.eType_TET) 
     2   io%e = "shapeFPrt only defined for tet elements"

      p => prt%get(ip)
      p%eID = 0
      iSb= prt%sb%id(p%x)
      b => prt%sb%box(iSb(1))


      do ii=1,size(b%els)

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
     2        msh%x(:,msh%IEN(a,b%els(ii)))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) +
     2        msh%x(:,msh%IEN(a,b%els(ii)))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + 
     2        msh%x(:,msh%IEN(a,b%els(ii)))*Nxi(3,a)
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
            prntx(a) = xiX(a,1)*(p%x(1) - msh%x(1,msh%IEN(4,b%els(ii))))
     2               + xiX(a,2)*(p%x(2) - msh%x(2,msh%IEN(4,b%els(ii))))
     3               + xiX(a,3)*(p%x(3) - msh%x(3,msh%IEN(4,b%els(ii))))
         END DO
      END IF

      ! Finding shape function values at particle coordinates
      N(1) = prntx(1)
      N(2) = prntx(2)
      N(3) = prntx(3)
      N(4) = 1 - prntx(1) - prntx(2) - prntx(3)

      ! Checking if all shape functions are positive
      IF (ALL(N.ge.0D0)) then
         p%eID=b%els(ii)
         EXIT
      END IF
         
      end do

      ! If it loops through everything and doesn't yield a positive shape function,
      ! the particle is outside the domain
      if (p%eID.eq.0) io%e = 'outside domain'

      p%eID = e
      p%N = N
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
           prt%dat = p
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
         e = prt%shapeF(ip,msh)
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
            fvel(ii) = fvel(ii) + u%A%v(ii,msh%IEN(jj,p%eID))*p%N(jj)
         end do
      end do

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
      SUBROUTINE collidePrt(prt,id1,id2,dtp)
      CLASS(prtType), INTENT(IN), TARGET :: prt
      INTEGER, INTENT(IN) :: id1, id2
      REAL(KIND=8), INTENT(IN) :: dtp
      TYPE(pRawType), POINTER :: p1,p2

!     Calculating distance coefficient
      REAL(KIND=8) :: a, b, c, d, e, f, qa, qb, qc, zeros(2), tcr
      REAL(KIND=8) :: n1(nsd), n2(nsd), t1(nsd), t2(nsd)
      REAL(KIND=8) :: vpar1, vpar2, vperp1, vperp2, dp, mp, rho, k
!     Coefficients to make calculating parallel/perp vel easier
      REAL(KIND=8) :: pa, pb

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
      p1%x = p1%u*tcr + p1%x
      p2%x = p2%u*tcr + p2%x

      ! Vector parallel and pependicular to collision tangent line
      n1 = (p1%x - p2%x)/((dp + dp)/2)
      n2 = -n1
      t1 = cross2(cross2(n1,p1%u),n1)
      t2 = cross2(cross2(n2,p2%u),n2)

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

      !!! Needs to be extended for multiple collisions per time step (will probably be here)
      ! Advance particle the rest of the time step at this velocity.
      p1%x = p1%x + p1%u*(dtp - tcr)
      p2%x = p2%x + p2%u*(dtp - tcr)

      p1%collided = .true.
      p2%collided = .true.

      p1%remdt = dtp-tcr
      p2%remdt = dtp-tcr

      END SUBROUTINE collidePrt

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
      SUBROUTINE solvePrt(prt, ns)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
      CLASS(eqType), INTENT(IN) :: ns
      INTEGER ip, a, e, eNoN, Ac,i ,j, k,l, subit
      REAL(KIND=8) rhoF,
     2   mu, mp,
     3   rhoP, N(nsd+1), Ntmp(nsd+1)
      TYPE(gVarType), POINTER :: u
      TYPE(sbType), POINTER :: sb
      TYPE(mshType), POINTER :: lM
      TYPE(pRawType), POINTER :: p(:)
      TYPE(pRawType), ALLOCATABLE :: tmpp
!     Particle/fluid Parameters
      REAL(KIND=8) :: g(nsd), rho, dtp,maxdtp,sbdt(nsd), dp, taup
!     Derived from flow
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd)
!     RK2 stuff
      REAL(KIND=8) :: prtxpred(nsd), pvelpred(nsd)

      ALLOCATE(p(prt%n))
      lM => ns%dmn%msh(1)
      u =>  ns%var(1)
      p =>  prt%dat
      sb => prt%sb
      rhoF  = ns%mat%rho()
      rhoP  = prt%mat%rho()
      mu    = ns%mat%mu()
      dp    = ns%mat%D()

      ! Particle relaxation time
      taup=rhop * dp**2D0/mu/18D0

            ! Gravity
!!!!! Find where grav actually is? (maybe mat%body forces)
      g=0D0

!     Initialize if haven't yet
      IF(.NOT.prt%crtd) THEN
            CALL prt%new(2)
      END IF

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

      DO l = 1,subit
      DO i = 1,prt%n
!        Find which searchbox particle is in
         p(i)%sbID = sb%id(p(i)%x)
!        Get shape functions/element of particle
         N = prt%shapeF(i, lM)
!        Get drag acceleration
         apd = prt%drag(i,ns)
!        Total acceleration (just drag and buoyancy now)
         apT = apd + g*(1D0 - rhoF/rhoP)
!        2nd order advance (Heun's Method)
!        Predictor
         pvelpred = p(i)%u + dtp*apT
         prtxpred = p(i)%x + dtp*p(i)%u
         tmpp%u   = pvelpred
         tmpp%x   = prtxpred
!        Find which searchbox prediction is in
         tmpp%sbID = sb%id(tmpp%x)
!        Get shape functions/element of prediction
         Ntmp = prt%shapeF(i, lM)         
!        Get drag acceleration of predicted particle
         apdpred = prt%drag(i,ns)
         apTpred = apdpred + g*(1D0 - rhoF/rhoP)
!        Corrector
         p(i)%u = p(i)%u + 0.5D0*dtp*(apT+apTpred)         
      END DO

      !!! REALLY need to figure out how to make it so it only checks for collisions in the same searchbox
      DO i=1,prt%n
!        Collisions
         DO j=1,prt%n
!        Check if the particle collides with any other particles and hasn't collided. Advance if so.
            IF ((i.ne.j).and.(.not.(p(i)%collided)))
     2          CALL prt%collide(i,j,dtp)
         ENDDO

!        If particles haven't collided, advance by vel*dtp
         IF (.not.(p(i)%collided)) THEN
            p(i)%x = dtp*p(i)%u + p(i)%x
         ELSE
            p(i)%collided = .false.
         END IF

         print *, p(i)%x,time
      END DO
      END DO


! mpiifort -O3 -module obj -I../memLS/include -I../cplBC/include -c src/PRT.f -o obj/PRT.o

      RETURN
      END SUBROUTINE solvePrt
      
      END MODULE PRTMOD
