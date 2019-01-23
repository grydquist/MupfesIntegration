      MODULE PRTMOD
      USE EQMOD
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: OLD=1, CUR=2, NEW=3, RHS=4, IMP=4

!     Size of the container is extended by this factor if it's
!     completely filled
      INTEGER, PARAMETER :: prtExtFac = 2

      TYPE pDisType
!        Number of particles to be seeded in the domain
         INTEGER :: n = 0
!        Particle release rate at inlet
         REAL(KIND=8) :: nd = 0D0
!        Uniform particle diameter
         REAL(KIND=8) dia
!        Uniform particle density
         REAL(KIND=8) rho
!        Restitution coefficient
         REAL(KIND=8) :: k = 1D0
!        Total mass of particles in this class
         REAL(KIND=8) m
!        The average velocity of particle in this class
         REAL(KIND=8) u(3)
!        Total kinetic energy 
         REAL(KIND=8) E
!        Collision frequency
         REAL(KIND=8) cf
!        RMS of relative velocity of colliding particles
         REAL(KIND=8) ur
!        Particle as tracers
         LOGICAL :: tracer = .FALSE.
!        Keep particles position fixed
         LOGICAL :: fixP = .FALSE.
      END TYPE pDisType


!!! may want to make another subtype later just called box (or use stack types) that
!!! contains individual box's + dimensions/elementss, then sb is just the collection

!     Search box type
      TYPE sbType
!        Created?
         LOGICAL :: crtd = .FALSE.
!     Number of boxes in each direction
         INTEGER n(3)
!     Size of boxes
         REAL(KIND=8) :: step(3)
!      Searchbox dimensions
         REAL(KIND=8), ALLOCATABLE :: dim(:,:)
!      Elements contained in searchbox
         INTEGER, ALLOCATABLE :: els(:,:)
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
!     Particle distribution ID
         INTEGER(KIND=8) :: dID=0
!     Position
         REAL(KIND=8) x(3,3)
!     Velocity
         REAL(KIND=8) u(3,4)
      END TYPE pRawType

      TYPE prtType
      !        Created?
         LOGICAL :: crtd = .FALSE.
!        Maximum capacity of container
         INTEGER :: maxN = 1024
!        Current size
         INTEGER :: n = 0
!        Number of particle distributions (classes)
         INTEGER :: nDis = 0
!        Search boxes to find element hosting particles
         TYPE(sbType) :: sb
!        Array pointing to full and empty entries
         INTEGER, ALLOCATABLE :: ptr(:)
!        Data
         TYPE(pRawType), ALLOCATABLE :: dat(:)
!        Particle distribution properties
         TYPE(pDisType), ALLOCATABLE :: dis(:)
!        Material properties
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
!     Advance all particles one time step
         PROCEDURE :: solve => solvePrt
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
      CLASS(sbType), INTENT(OUT):: sb
      INTEGER :: ii,jj,cnt2,kk
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), elbox(2*nsd,msh%nEl)
      INTEGER, ALLOCATABLE :: sbel(:)

 !!     ! Will need to replace eventually with something that targets element size
      INTEGER :: split(nsd)
      split=(/10,10,10/)
      sb%n=split

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sb%dim(2*nsd,sb%n(1)*sb%n(2)*sb%n(3)))

      ALLOCATE(sbel(msh%nEl))
      ALLOCATE(sb%els(sb%n(1)*sb%n(2)*sb%n(3),msh%nEl))

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
         sb%dim(1,(seq1+ii-1)) = MINVAL(msh%x(1,:)) 
     2       + sb%step(1)*(ii-1)/2
      end do

      ! Direction 2
      do ii=1,sb%n(2)
         sb%dim(3,seq2+(ii-1)*split(1)) = MINVAL(msh%x(2,:))
     2       + sb%step(2)*(ii-1)/2
      end do

      ! Direction 3
      do ii=1,sb%n(3)
         sb%dim(5,seq3+(ii-1)*split(1)*split(2))=MINVAL(msh%x(3,:))
     2        +sb%step(3)*(ii-1)/2
      end do

      sb%dim(2,:) = sb%dim(1,:) + sb%step(1)
      sb%dim(4,:) = sb%dim(3,:) + sb%step(2)
      sb%dim(6,:) = sb%dim(5,:) + sb%step(3)

      ! Making boxes surrounding elements
      do ii=1,msh%Nel
         do jj=1,nsd
            elbox(2*jj-1,ii) = MINVAL(msh%x(jj,msh%IEN(:,ii)))
            elbox(2*jj  ,ii) = MAXVAL(msh%x(jj,msh%IEN(:,ii)))
         end do
      end do

      do ii=1,split(1)*split(2)*split(3)
         cnt2=1
         sbel=0
         do jj=1,msh%Nel
            ! Check if elements are completely outside searchbox
            do kk=1,nsd
               ! Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).lt.sb%dim(ii,2*kk  )) cycle
               if (elbox(2*kk  ,jj).gt.sb%dim(ii,2*kk-1)) cycle
            end do

            sbel(cnt2) = jj
            cnt2=cnt2+1
            !end if
         end do
         sb%els(ii,:)=sbel
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
      xzero(1) = x(1) - minval(sb%dim(1,:))
      xzero(2) = x(2) - minval(sb%dim(3,:))
      xzero(3) = x(3) - minval(sb%dim(5,:))

      ! Find which searchbox the particle is in
      ! Number of searchbox steps in x,y,and z
      xsteps=FLOOR(xzero/sb%step)
      ! furthest searchbox in front
      iSb(1) = xsteps(1) + sb%n(1)*xsteps(2)+
     2   sb%n(1)*sb%n(2)*xsteps(3)+1
      ! previous sb in x
      iSb(2) = iSb(1)-1
      ! previous sb's in y
      iSb(3) = iSb(1)-sb%n(1)
      iSb(4) = iSb(3)-1
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
!     First checks to see if the particle is still in this cell,
!     otherwise will search all the elements in the search box. Returns
!     -1 if particle is outisde of the mesh

!!! Doesn't work
      FUNCTION shapeFPrt(prt, ip, N,msh) RESULT(e)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN) :: prt
      INTEGER, INTENT(IN) :: ip
      REAL(KIND=8), INTENT(OUT) :: N(4)
      TYPE(mshType), INTENT(IN) :: msh
      INTEGER e

      INTEGER, PARAMETER :: TIME_POINT=CUR, iM=1

      LOGICAL NOX
      INTEGER ie , iSb(2**nsd)
      TYPE(pRawType), POINTER :: p

      IF (msh%eType.NE.eType_TET) 
     2   io%e = "shapeFPrt only defined for tet elements"

      p => prt%get(ip)
      IF (p%eID .GT. 0) THEN
         IF (NOX(msh,p%x(:,TIME_POINT),p%eID,N)) RETURN
      END IF

      iSb = prt%sb%id(p%x(:,TIME_POINT))
!      DO ie=1, prt%sb%c(iSb)%n
!         e = prt%sb%c(iSb)%v(ie)
!         IF (NOX(msh,p%x(:,TIME_POINT),e,N)) EXIT
!      END DO
      e=1
      !IF (ie .GT. prt%sb%c(iSb)%n) e = -1
      p%eID = e

      RETURN
      END FUNCTION shapeFPrt
!--------------------------------------------------------------------
      SUBROUTINE seedPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
      
      INTEGER ip, iDis
      TYPE(pDisType), POINTER :: pd
      TYPE(pRawType) p

      CALL RSEED(cm%id())
      DO iDis=1, prt%nDis
         pd => prt%dis(iDis)
         DO ip=1, pd%n
            p%u    = 0D0
            p%dID  = iDis
            p%pID  = INT(ip,8)
            CALL RANDOM_NUMBER(p%x(:,OLD))
            IF (nsd .EQ. 2) p%x(3,OLD) = 0D0
            p%x(:,OLD) = (p%x(:,OLD) - (/0.5D0,0.5D0,0D0/))*2D0
            p%x(:,CUR) = p%x(:,OLD)
            p%x(:,NEW) = p%x(:,OLD)
            CALL prt%add(p)
         END DO
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
      prt%nDis = 1
      ALLOCATE(prt%dis(prt%nDis))
      prt%dis(1)%n = 5
      CALL prt%seed()
      CALL prt%sb%new(msh)



      DO ip=1, prt%n
         p => prt%get(ip)
         e = prt%shapeF(ip,N,msh)
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
      SUBROUTINE solvePrt(prt, ns)!, Rm)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
      CLASS(eqType), INTENT(IN) :: ns
      ! TYPE(vecType), INTENT(OUT) :: Rm

      INTEGER ip, a, e, eNoN, Ac
      REAL(KIND=8) ug(3), ap(3), f(3), up(3), um, rhoF,
     2   mu, fL(3), mp, fD(3), Cd, N(4), g(3),
     3   fT(3), us(3), rhoP
      TYPE(pDisType), POINTER :: pd
      TYPE(pRawType), POINTER :: p
      TYPE(gVarType), POINTER :: u
      TYPE(mshType), POINTER :: lM

      u => ns%var(1)
      !u%A%v(i,Ac) ! velocity of node Ac in direction i
      lM => ns%dmn%msh(1)
      !DO g=1, lM%nG
!            IF (g.EQ.1.OR..NOT.lM%lShpF) CALL lM%dNdx(g, xl, Nx, J, ks)
!            IF (ISZERO(J)) io%e = "Jac < 0 @ element "//e
      !      w = lM%w(g)*J
      !      N = lM%N(:,g)
      !END DO
      rhoF  = ns%mat%rho()
      rhoP  = prt%mat%rho()
      mu    = ns%mat%mu()


      !rhsRU = 0D0
      DO ip=1, prt%n

!     Particle mass
         mp = pi/6D0*rhoP*pd%dia**3D0
!     Particle velocity
         up = p%u(:,CUR) + p%u(:,IMP)
!     Actual slip velocity
         us = up - ug
!     Relative velocity magnitude
         um = SQRT(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
!     Finite Re effect
         Cd = 1D0 + 0.15D0*(um*pd%dia*rhoF/mu)**0.687D0
!     Drag force
         fD = -3D0*pi*mu*pd%dia*us*Cd
!     Particle acceleration
         ap = fT/mp
!     Total acceleration
         ap = ap + g*(1D0 - rhoF/rhoP)
!     No acceleration if particle is fixed
         IF (pd%fixP) ap = 0D0
!     Tracers have same velocity as fluid
         IF (pd%tracer) THEN
            ap         = 0D0
            p%u(:,CUR) = ug
            p%u(:,NEW) = ug
         END IF

!     Resetting impact function. To be set by the collisons calls
         p%u(:,IMP) = 0D0
      END DO
!     Particle-particle collisions
c      IF (prt%fwc) CALL prt%collision()
!     And this is particle wall collision
c      CALL prt%wallCorrect()
!     Doing communication and updating indices associated with particle
!     locations on the mesh
c      CALL prt%comu()
c      CALL prt%upInd()

      RETURN
      END SUBROUTINE solvePrt
      
      END MODULE PRTMOD


! Need to fix whole prttype, then integrate all the shit into prtsolve
! prtsolve will likely use all other subroutines in readvtk