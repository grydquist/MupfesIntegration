      MODULE PRTMOD

      USE COMMOD
      USE ALLFUN
      USE CHNLMOD

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

!     Search box type
      TYPE sbType
!        Created?
         LOGICAL :: crtd = .FALSE.
!     Number of boxes in each direction
         INTEGER n(3)
!     The lower bound of searched domain in each direction
         REAL(KIND=8) x(3)
!     The size of search box in each direction
         REAL(KIND=8) d(3)
!     Container to keep track of elements within each box
         TYPE(stackType), ALLOCATABLE :: c(:)
      CONTAINS
!     Sets up the search boxes pertaining to msh(1)
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
      
      IF (prt%sb%crtd) THEN
         DO iSb=1, prt%sb%n(1)*prt%sb%n(2)*prt%sb%n(3)
            CALL prt%sb%c(iSb)%free()
         END DO
         DEALLOCATE(prt%sb%c)
      END IF

      RETURN
      END SUBROUTINE freePrt
!--------------------------------------------------------------------
      PURE SUBROUTINE newSb(sb)
      IMPLICIT NONE
      CLASS(sbType), INTENT(INOUT) :: sb

!     Assuming everything is located in mesh 1
      INTEGER, PARAMETER :: iM=1
!     Element to search box size ratio
      REAL(KIND=8), PARAMETER :: R=1D0

      INTEGER iSb, e, i, j, k, iL(nsd), iH(nsd)
      REAL(KIND=8) xmax(nsd), xmin(nsd), l(nsd)
      
      IF (sb%crtd) RETURN

!     Maximum element position in each direction
      xmax = MAXVAL(x,2) 
      xmax = xmax + ABS(xmax)*eps + eps
!     Minimum element position in each direction (lower bound of SB)
      sb%x = MINVAL(x,2)
      sb%x = sb%x - ABS(sb%x)*eps - eps
!     Finding the size of search box such the cost of computation is
!     minimized. So we target 1 elements per search box
      l = xmax - sb%x
      IF (nsd .EQ. 2) THEN
!     If we want 1 particle per box, what should be the size of the box
         sb%d = SQRT(l(1)*l(2)/msh(iM)%gnEl)/R
!     This corresponds to an estimate of number of boxes
         sb%n = MAX(INT(l/sb%d),(/1,1/))
         ALLOCATE(sb%c(sb%n(1)*sb%n(2)))
      ELSE
         sb%d = ((l(1)*l(2)*l(3)/msh(iM)%gnEl)**(1D0/3D0))/R
!     This corresponds to an estimate of number of boxes
         sb%n = MAX(INT(l/sb%d),(/1,1,1/))
         ALLOCATE(sb%c(sb%n(1)*sb%n(2)*sb%n(3)))
      END IF
!     So this is the final box size
      sb%d = l/sb%n

!     Populating boxes
      DO e=1, msh(iM)%nEl
         DO i=1, nsd
            xmax(i) = MAXVAL(x(i,msh(iM)%IEN(:,e)))
            xmin(i) = MINVAL(x(i,msh(iM)%IEN(:,e)))
         END DO
         iL = CEILING((xmin - sb%x)/sb%d)
         iH = CEILING((xmax - sb%x)/sb%d)
         DO i=iL(1), iH(1)
            DO j=iL(2), iH(2)
               IF (nsd .EQ. 2) THEN
                  iSb = (j-1)*sb%n(1) + i
                  CALL sb%c(iSb)%push(e)
               ELSE
                  DO k=iL(3), iH(3)
                     iSb = ((k-1)*sb%n(2) + j-1)*sb%n(1) + i
                     CALL sb%c(iSb)%push(e)
                  END DO
               END IF
            END DO
         END DO
      END DO
      sb%crtd = .TRUE.

      RETURN
      END SUBROUTINE newSb
!--------------------------------------------------------------------
!     Returns the ID of a searchbox that contains point x
      PURE FUNCTION idSB(sb,x) RESULT(iSb)
      IMPLICIT NONE
      CLASS(sbType), INTENT(IN) :: sb
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER iSb

      INTEGER i(3)
      
      i = CEILING((x-sb%x)/sb%d)
      IF (nsd .EQ. 2) THEN
         iSb = (i(2) - 1)*sb%n(1) + i(1)
      ELSE
         iSb = ((i(3) - 1)*sb%n(2) + i(2) - 1)*sb%n(1) + i(1)
      END IF
      
      RETURN
      END FUNCTION idSB
!--------------------------------------------------------------------
!     First checks to see if the particle is still in this cell,
!     otherwise will search all the elements in the search box. Returns
!     -1 if particle is outisde of the mesh
      FUNCTION shapeFPrt(prt, ip, N) RESULT(e)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN) :: prt
      INTEGER, INTENT(IN) :: ip
      REAL(KIND=8), INTENT(OUT) :: N(4)
      INTEGER e

      INTEGER, PARAMETER :: TIME_POINT=CUR, iM=1

      LOGICAL NOX
      INTEGER ie , iSb
      TYPE(pRawType), POINTER :: p

      IF (msh(iM)%eType.NE.eType_TET) 
     2   err = "shapeFPrt only defined for tet elements"

      p => prt%get(ip)
      IF (p%eID .GT. 0) THEN
         IF (NOX(msh(iM),p%x(:,TIME_POINT),p%eID,N)) RETURN
      END IF

      iSb = prt%sb%id(p%x(:,TIME_POINT))
      DO ie=1, prt%sb%c(iSb)%n
         e = prt%sb%c(iSb)%v(ie)
         IF (NOX(msh(iM),p%x(:,TIME_POINT),e,N)) EXIT
      END DO
      IF (ie .GT. prt%sb%c(iSb)%n) e = -1
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
 
      eNoN = msh(1)%eNoN
!     Initializing particles
      CALL prt%new()
      prt%nDis = 1
      ALLOCATE(prt%dis(prt%nDis))
      prt%dis(1)%n = 5
      CALL prt%seed()
      CALL prt%sb%new()



      DO ip=1, prt%n
         p => prt%get(ip)
         e = prt%shapeF(ip,N)
         up = 0D0
         DO a=1, eNoN
            Ac = msh(1)%IEN(a,e)
            up = up + Yn(1:3,Ac)*N(a)
         END DO
      END DO


      call prt%free()

      RETURN
      END SUBROUTINE PRTTRANS
!--------------------------------------------------------------------
      SUBROUTINE corXUPrt(prt)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
c      TYPE(vecType), INTENT(INOUT) :: rhsRU

      INTEGER ip, a, e, eNoN, Ac
      REAL(KIND=8) ug(3), ap(3), f(3), up(3), um, rho,
     2   mu, fL(3), mp, fD(3), Cd, N(4), g(3),
     3   fT(3), us(3)
      TYPE(pDisType), POINTER :: pd
      TYPE(pRawType), POINTER :: p

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      mu   = eq(cEq)%dmn(cDmn)%prop(viscosity)
      g(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      g(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      g(3) = eq(cEq)%dmn(cDmn)%prop(f_z) 

      !rhsRU = 0D0
      DO ip=1, prt%n
         p  => prt%get(ip)
         pd => prt%dis(p%dID)
         
         e = prt%shapeF(ip,N)
         ug = 0D0
         DO a=1, eNoN
            Ac = msh(1)%IEN(a,e)
            ug = ug + Yn(1:3,Ac)*N(a)
         END DO

!     Particle mass
         mp = pi/6D0*pd%rho*pd%dia**3D0
!     Particle velocity
         up = p%u(:,CUR) + p%u(:,IMP)
!     Actual slip velocity
         us = up - ug
!     Relative velocity magnitude
         um = SQRT(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
!     Finite Re effect
         Cd = 1D0 + 0.15D0*(um*pd%dia*rho/mu)**0.687D0
!     Drag force
         fD = -3D0*pi*mu*pd%dia*us*Cd
!     Particle acceleration
         ap = fT/mp
!     Total acceleration
         ap = ap + g*(1D0 - rho/pd%rho)
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
      END SUBROUTINE corXUPrt
      
      END MODULE PRTMOD

!     for checking particle positions before/after calling NAtx
      isInside = .FALSE.
      DO i=1, nsd
         IF (x(i) .GT. MAXVAL(xl(i,:))) RETURN
         IF (x(i) .LT. MINVAL(xl(i,:))) RETURN
      END DO
 
      isInside = ALL(N .GT. 0D0)


