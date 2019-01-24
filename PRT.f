      MODULE PRTMOD
      USE EQMOD
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: OLD=1, CUR=2, NEW=3, RHS=4, IMP=4

!     Size of the container is extended by this factor if it's
!     completely filled
      INTEGER, PARAMETER :: prtExtFac = 2

!!! may want to make another subtype later just called box (or use stack types) that
!!! contains individual box's + dimensions/elementss, then sb is just the collection

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
!     Position
         REAL(KIND=8) x(3)
!     Velocity
         REAL(KIND=8) u(3)
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
         do jj=1,msh%Nel
            ! Check if elements are completely outside searchbox
            do kk=1,nsd
               ! Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).lt.sb%box(ii)%dim(2*kk  )) cycle
               if (elbox(2*kk  ,jj).gt.sb%box(ii)%dim(2*kk-1)) cycle
            end do

            sbel(cnt2) = jj
            cnt2=cnt2+1
            !end if
         end do
         ALLOCATE(sb%box(ii)%dim(2*nsd))
         ALLOCATE(sb%box(ii)%els(msh%nEl))
         sb%box(ii)%els=sbel
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
!     otherwise will search all the elements in the search box.

      FUNCTION shapeFPrt(prt, ip, N,msh) RESULT(e)
      IMPLICIT NONE
      CLASS(prtType), INTENT(IN),TARGET :: prt
      INTEGER, INTENT(IN) :: ip
      REAL(KIND=8), INTENT(OUT) :: N(4)
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
      p%eID=0
      iSb = prt%sb%id(p%x)
      b => prt%sb%box(iSb(1))


      do ii=1,msh%nEl

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
      N(4) = 1 - prntx(1) - prntx(2)
     2 - prntx(3)

      ! Checking if all shape functions are positive
      IF (ALL(N.gt.0D0)) then
         p%eID=b%els(ii)
         EXIT
      END IF
         
      end do

      ! If it loops through everything and doesn't yield a positive shape function,
      ! the particle is outside the domain
      if (p%eID.eq.0) print *, 'outside domain'

      p%eID = e
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
      prt%nDis = 1
      !ALLOCATE(prt%dis(prt%nDis))
      !prt%dis(1)%n = 5
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
      SUBROUTINE solvePrt(prt, ns)
      IMPLICIT NONE
      CLASS(prtType), INTENT(INOUT), TARGET :: prt
      CLASS(eqType), INTENT(IN) :: ns

      INTEGER ip, a, e, eNoN, Ac
      REAL(KIND=8) ug(3), ap(3), f(3), up(3), um, rhoF,
     2   mu, fL(3), mp, fD(3), Cd, N(4), g(3),
     3   fT(3), us(3), rhoP
      TYPE(gVarType), POINTER :: u
      TYPE(mshType), POINTER :: lM

      lM => ns%dmn%msh(1)
      u => ns%var(1)

      IF(.NOT.prt%crtd) THEN
            CALL prt%new(2)
      END IF

      IF (.NOT.(prt%sb%crtd)) THEN
         CALL prt%sb%new(lM)
      END IF



      !u%A%v(i,Ac) ! velocity of node Ac in direction i

      !DO g=1, lM%nG
!            IF (g.EQ.1.OR..NOT.lM%lShpF) CALL lM%dNdx(g, xl, Nx, J, ks)
!            IF (ISZERO(J)) io%e = "Jac < 0 @ element "//e
      !      w = lM%w(g)*J
      !      N = lM%N(:,g)
      !END DO
      rhoF  = ns%mat%rho()
      rhoP  = prt%mat%rho()
      mu    = ns%mat%mu()

      RETURN
      END SUBROUTINE solvePrt
      
      END MODULE PRTMOD



! prtsolve will likely use all other subroutines in readvtk