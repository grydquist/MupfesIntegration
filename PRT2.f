      MODULE PRT2MOD
      
      USE COMMOD
      USE ALLFUN
      USE CHNLMOD

      IMPLICIT NONE


      TYPE prt
      ! Properties
         REAL(KIND=8) :: mp, dp = 1D0, rhop = 0.1D0
      ! Flow characteristics
         REAL(KIND=8) :: x(nsd), vel(nsd),prntx(nsd),shps(eNoN), remdtp
      ! Searechbox/ element location
         INTEGER :: sbid(2**nsd), elid, near(Np)
      ! Collisions with other particles
         LOGICAL :: collided=.false.
      ! Restitution Coefficient
         REAL(KIND=8) :: k
      END TYPE prt


      ! Collection of particles
      type(prt) :: prts(Np)

      TYPE sb
      ! Size of searchbox
         REAL(KIND=8) :: step(nsd)
      ! Searchbox dimensions
         REAL(KIND=8) :: dim(nsd*2)
      ! Elements contained in searchbox
         INTEGER, ALLOCATABLE :: els(:)
      END TYPE sb

      ! Domain spliited into sb's
      type(sb), ALLOCATABLE :: sbdom(:)

      REAL, PARAMETER :: pi=3.1415926535897932384626433

      prts%mp = pi*prts%rhop/6D0*prts%dp**3D0

      CONTAINS

!#################################################################### SBDOMAIN
      SUBROUTINE SBDomain(x,split,sbdom)
      REAL(KIND=8), INTENT(IN) :: x(nsd,nNo)
      INTEGER, INTENT(IN) :: split(nsd)
      TYPE(sb), INTENT(OUT), ALLOCATABLE :: sbdom(:)
      INTEGER :: ii,jj,cnt2,kk
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), elbox(2*nsd,msh(1)%gnEl)
      INTEGER, ALLOCATABLE :: sbel(:)

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sbdom(split(1)*split(2)*split(3)))

      ALLOCATE(sbel(msh(1)%gnEl))

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(split(3)*split(2)),seq2(split(3)*split(1)),seq3(split(2)*split(1)))

      ! Domain ranges
      diff(1)=MAXVAL(x(1,:))-MINVAL(x(1,:))
      diff(2)=MAXVAL(x(2,:))-MINVAL(x(2,:))
      diff(3)=MAXVAL(x(3,:))-MINVAL(x(3,:))
      ! Size of sb
      do ii = 1,(split(2)*split(2)*split(3))
         sbdom(ii)%step=diff/((split+1)/2)
      end do

      seq1=(/(ii, ii=0, split(2)*split(3)-1, 1)/)*split(1)+1
      cnt2=0
      do ii=1,split(1)*split(3)
            seq2(ii)=ii+cnt2*(split(2)-1)*split(1)
            if (MOD(ii,split(1)).eq.0) cnt2=cnt2+1
      end do
      seq3=(/(ii, ii=0, split(1)*split(2)-1, 1)/)+1

      ! Allocating sbdim, such that they overlap by 50%
      do ii=1,split(1)
         sbdom(seq1+ii-1)%dim(1) = MINVAL(x(1,:)) + sbdom(1)%step(1)*(ii-1)/2
      end do

      do ii=1,split(2)
         sbdom(seq2+(ii-1)*split(1))%dim(3) = MINVAL(x(2,:)) + sbdom(1)%step(2)*(ii-1)/2
      end do

      do ii=1,split(3)
         sbdom(seq3+(ii-1)*split(1)*split(2))%dim(5)=MINVAL(x(3,:))+sbdom(1)%step(3)*(ii-1)/2
      end do

      sbdom%dim(2) = sbdom%dim(1) + sbdom(1)%step(1)
      sbdom%dim(4) = sbdom%dim(3) + sbdom(1)%step(2)
      sbdom%dim(6) = sbdom%dim(5) + sbdom(1)%step(3)

      ! Making boxes surrounding elements
      do ii=1,msh(1)%gNel
         do jj=1,nsd
            elbox(2*jj-1,ii) = MINVAL(x(jj,msh(1)%IEN(:,ii)))
            elbox(2*jj  ,ii) = MAXVAL(x(jj,msh(1)%IEN(:,ii)))
         end do
      end do

      do ii=1,split(1)*split(2)*split(3)
         cnt2=1
         sbel=0
         do jj=1,msh(1)%gNel
            ! Check if elements are completely outside searchbox
            do kk=1,nsd
               ! Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).lt.sbdom(ii)%dim(2*kk  )) cycle
               if (elbox(2*kk  ,jj).gt.sbdom(ii)%dim(2*kk-1)) cycle
            end do

            sbel(cnt2) = jj
            cnt2=cnt2+1
            !end if
         end do
         sbdom(ii)%els=sbel
      end do
      
      END SUBROUTINE SBDomain

!#################################################################### XSB
      ! Find all Searchboxes x is in
      SUBROUTINE xSB(myprt,sbdom,split)
      IMPLICIT NONE
      TYPE(sb), INTENT(IN), ALLOCATABLE :: sbdom(:)
      TYPE(prt), INTENT(INOUT) :: myprt
      INTEGER, INTENT(IN) :: split(nsd)
      REAL(KIND=8) :: xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Set domain back to zero
      xzero(1) = myprt%x(1) - minval(sbdom%dim(1))
      xzero(2) = myprt%x(2) - minval(sbdom%dim(3))
      xzero(3) = myprt%x(3) - minval(sbdom%dim(5))

      ! Find which searchbox the particle is in
      ! Number of searchbox steps in x,y,and z
      xsteps=FLOOR(xzero/sbdom(1)%step)
      ! furthest searchbox in front
      myprt%sbid(1) = xsteps(1)+split(1)*xsteps(2)+split(1)*split(2)*xsteps(3)+1
      ! previous sb in x
      myprt%sbid(2) = myprt%sbid(1)-1
      ! previous sb's in y
      myprt%sbid(3) = myprt%sbid(1)-split(1)
      myprt%sbid(4) = myprt%sbid(3)-1
      ! Next sb's in z (if available)
      if (nsd.eq.3) then
         myprt%sbid(5) = myprt%sbid(1) - split(1)*split(2)
         myprt%sbid(6) = myprt%sbid(5) - 1
         myprt%sbid(7) = myprt%sbid(5) - split(1)
         myprt%sbid(8) = myprt%sbid(7) - 1
      end if
      

      END SUBROUTINE xSB

!#################################################################### XEL
      ! Finds element particle of position x is in is in
      SUBROUTINE xEl(sbdom,myprt,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: x(nsd,nNo)
      TYPE(sb), INTENT(IN) :: sbdom
      TYPE(prt), INTENT(INOUT) :: myprt
      INTEGER :: ii,cnt,a
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,eNoN)
      cnt=1
      myprt%elid=0

      do ii=1,msh(1)%gnEl

      xXi = 0D0

      IF (nsd .EQ. 2) THEN
      !
      ! 2D not done
      !
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         
         DO a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         END DO


      ELSE

         DO a=1, eNoN
           xXi(:,1) = xXi(:,1) + x(:,msh(1)%IEN(a,sbdom%els(ii)))*Nxi(1,a)
           xXi(:,2) = xXi(:,2) + x(:,msh(1)%IEN(a,sbdom%els(ii)))*Nxi(2,a)
           xXi(:,3) = xXi(:,3) + x(:,msh(1)%IEN(a,sbdom%els(ii)))*Nxi(3,a)
         END DO
         
         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)&
     &       + xXi(1,2)*xXi(2,3)*xXi(3,1)&
     &       + xXi(1,3)*xXi(2,1)*xXi(3,2)&
     &       - xXi(1,1)*xXi(2,3)*xXi(3,2)&
     &       - xXi(1,2)*xXi(2,1)*xXi(3,3)&
     &       - xXi(1,3)*xXi(2,2)*xXi(3,1)

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
      myprt%prntx(a) = xiX(a,1)*(myprt%x(1) - x(1,msh(1)%IEN(4,sbdom%els(ii)))) + &
     &                 xiX(a,2)*(myprt%x(2) - x(2,msh(1)%IEN(4,sbdom%els(ii)))) + &
     &                 xiX(a,3)*(myprt%x(3) - x(3,msh(1)%IEN(4,sbdom%els(ii))))
         END DO
      END IF

      myprt%shps(1) = myprt%prntx(1)
      myprt%shps(2) = myprt%prntx(2)
      myprt%shps(3) = myprt%prntx(3)
      myprt%shps(4) = 1 - myprt%prntx(1) - myprt%prntx(2) - myprt%prntx(3)

      IF (ALL(myprt%shps.gt.0D0)) then
         myprt%elid=sbdom%els(ii)
         EXIT
      END IF
         
      end do

      if (myprt%elid.eq.0) print *, 'outside domain'

      END SUBROUTINE xEl

!#################################################################### PRTDRAG
      ! Find acceleration on one particle from drag
      SUBROUTINE prtDrag(myprt,apd,taupo)
      TYPE(prt), INTENT(IN) :: myprt
      REAL(KIND=8), INTENT(OUT) :: apd(nsd)
      REAL(KIND=8) :: fvel(nsd),taup
      INTEGER ii,jj
      REAL(KIND=8), INTENT(OUT), OPTIONAL :: taupo

      !fluid Parameters
      REAL(KIND=8) :: rho, mu
      ! Derived from flow
      REAL(KIND=8) :: fSN, magud, Rep, relvel(nsd)


      ! Fluid density
      rho=1D0
      ! Fluid viscosity
      mu=0.01D0

      ! Particle relaxation time
      taup=myprt%rhop * myprt%dp**2D0/mu/18D0
      if (present(taupo)) taupo=taup

      ! Interpolate velocity at particle point
      fvel=0D0
      do ii=1,nsd
         do jj=1,eNoN
            fvel(ii) = fvel(ii) + Yn(ii,msh(1)%IEN(jj,myprt%elid))*myprt%shps(jj)
         end do
      end do

      ! Relative velocity
      relvel = fvel-myprt%vel
      ! Relative velocity magnitude
      magud = SUM(relvel**2D0)**0.5D0
      ! Reynolds Number
      Rep = myprt%dp*magud*myprt%rhop/mu
      ! Schiller-Neumann (finite Re) correction
      fSN = 1D0 + 0.15D0*Rep**0.687D0
      ! Stokes corrected drag force
      apD = fSN/taup*relvel
      
      END SUBROUTINE prtDrag

!#################################################################### PRTADVANCE
      ! Advance 1 Particle through flow
      SUBROUTINE prtAdvance(myprt)
      IMPLICIT NONE
      TYPE(prt), INTENT(INOUT) :: myprt
      TYPE(prt) :: tmpprt
      INTEGER ii

      !Particle/fluid Parameters
      REAL(KIND=8) :: g(nsd), rho, dtp,maxdtp,sbdt(nsd)
      ! Derived from flow
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd), taup
      ! RK2 stuff
      REAL(KIND=8) :: prtxpred(nsd), pvelpred(nsd)

      tmpprt = myprt

      ! Gravity
      g=0D0
      !g(3)=1D0
      ! Fluid density
      rho=1D0

      !! Time step (still need to match to overall flow solver)
      ! Maxdt for overall solver
      maxdtp = 0.01D0

      ! Get drag acceleration/particle relaxation time
      CALL prtDrag(myprt,apd,taup)

      ! Separate into sb sizes (could be optimized here)
      do ii=1,nsd
         sbdt(ii)=(sbdom(1)%dim(2*ii)-sbdom(1)%dim(2*ii-1))/abs(myprt%vel(ii))
      end do
      
      ! dtp is minimum between time to travel half searchbox, 1/10 relaxation time, and maxdtp
      dtp = min(maxdtp,0.5*minval(sbdt),taup/10)

      ! Total acceleration (just drag and buoyancy now)

      apT = apd + g*(1D0 - rho/myprt%rhop)

      ! 2nd order advance (Heun's Method)
      ! Predictor
      pvelpred = myprt%vel + dtp*apT
      prtxpred = myprt%x + dtp*myprt%vel

      tmpprt%vel = pvelpred
      tmpprt%x   = prtxpred

      CALL xSB(tmpprt,sbdom,split)
      CALL xEl(sbdom(tmpprt%sbid(1)),tmpprt,x)
      CALL prtDrag(tmpprt,apdpred)

      apTpred = apdpred + g*(1D0 - rho/myprt%rhop)

      ! Corrector
      myprt%vel = myprt%vel + 0.5D0*dtp*(apT+apTpred)

      ! Collisions work. Just need to figure out how to change velocity of both particles
      !! Maybe get only velocities for particles, bring those out, then advance all particles through flow in
      !!    collision solver?
      
      !prtx = prtx + 0.5D0*dtp*(pvel+pvelpred)
      time=time+dtp

      END SUBROUTINE prtAdvance

!#################################################################### PRTCOLLIDE
      ! Detects and enacts collisions
      !! Only between particles right now
      SUBROUTINE prtCollide(prt1,prt2,dtp)
      IMPLICIT NONE
      TYPE(prt), INTENT(INOUT) :: prt1, prt2
      REAL(KIND=8), INTENT(IN) :: dtp

      ! Calculating distance coefficient
      REAL(KIND=8) :: a, b, c, d, e, f, qa, qb, qc, zeros(2), tcr
      REAL(KIND=8) :: n1(nsd), n2(nsd), t1(nsd), t2(nsd)
      REAL(KIND=8) :: vpar1, vpar2, vperp1, vperp2
      ! Coefficients to make calculating parallel/perp vel easier
      REAL(KIND=8) :: pa, pb

      prt1%collided = .false.
      prt1%collided = .false.

      ! First, check if particles will collide at current trajectory
      a = prt1%x(1)   - prt2%x(1)
      b = prt1%vel(1) - prt2%vel(1)
      c = prt1%x(2)   - prt2%x(2)
      d = prt1%vel(2) - prt2%vel(2)
      if(nsd.eq.3) then
         e = prt1%x(3)   - prt2%x(3)
         f = prt1%vel(3) - prt2%vel(3)
      else
         e=0D0
         f=0D0
      end if
      
      qa = b**2D0 + d**2D0 + f**2D0
      qb = 2D0*(a*b + c*d +e*f)
      qc = a**2D0 + c**2D0 + e**2D0 - ((prt1%dp + prt2%dp)/2D0)**2D0

      ! Imaginary zeros means particles won't collide
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
      prt1%x = prt1%vel*tcr + prt1%x
      prt2%x = prt2%vel*tcr + prt2%x

      ! Vector parallel and pependicular to collision tangent line
      n1 = (prt1%x - prt2%x)/((prt1%dp + prt2%dp)/2)
      n2 = -n1
      t1 = cross(cross(n1,prt1%vel),n1)
      t2 = cross(cross(n2,prt2%vel),n2)

      ! Rare case with no perpendicular velocity
      if (ANY(ISNAN(t1))) t1 = 0D0
      if (ANY(ISNAN(t2))) t2 = 0D0
      
      ! Get precollision parallel and perpendicular velocities
      vperp1 = sum(t1*prt1%vel)
      vpar1  = sum(n1*prt1%vel)
      vperp2 = sum(t2*prt2%vel)
      vpar2  = sum(n2*prt2%vel)

      ! Note that perpendicular velocities don't change, so we only need to calculate parallel
      pa = prt1%mp*vpar1 - prt2%mp*vpar2
      pb = (-vpar1 - vpar2)*prt1%k

      vpar2 = (pa - prt1%mp*pb)/(prt1%mp + prt2%mp)
      vpar1 = pb + vpar2
      vpar2 = -vpar2

      ! V here is split into just two velocities, so just add them as vector

      prt1%vel = vpar1*n1 + vperp1*t1
      prt2%vel = vpar2*n2 + vperp2*t2

      !!! Needs to be extended for multiple collisions per time step (will probably be here)
      !! Doesn't work if particle is still because of cross products
      ! Advance particle the rest of the time step at this velocity.
      prt1%x = prt1%x + prt1%vel*(dtp - tcr)
      prt2%x = prt2%x + prt2%vel*(dtp - tcr)

      prt1%collided = .true.
      prt2%collided = .true.

      prt1%remdtp = dtp-tcr
      prt2%remdtp = dtp-tcr


      END SUBROUTINE prtCollide

!#################################################################### CROSS
      ! I use cross products a couple times above. Also normalizes to unit vector
      FUNCTION cross(v1,v2)
      IMPLICIT NONE
      REAL(KIND=8) :: cross(nsd)
      REAL(KIND=8), INTENT(IN) :: v1(nsd), v2(nsd)

      cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
      cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
      cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

      cross = cross/sqrt(cross(1)**2D0 + cross(2)**2D0 + cross(3)**2D0)
      
      END FUNCTION cross




!!!! right now I need to get x, vel, nEl, nNo, dtp
! nEl = msh(1)% gNel
! x is just
! I only use Nno for declaring, likely won't need to declare x/vel anymore
! IEN = msh(1)%IEN
! vel... Yn perhaps?

!! I removed the vel and x arguments from prtdrag and down

!! Right now I just need to see about compilation

! Here's how you handle collisions:
!           do i=1,Np
!               ! Collisions
!            do j=1,Np
!               ! Check if the particle collides with any other particles and hasn't collided. Advance if so.
!               if ((i.ne.j).and.(.not.(prts(i)%collided))) CALL prtCollide(prts(i),prts(j),0.01D0)
!            end do
!
!            ! If particles haven't collided, advance by vel*dtp
!            if (.not.(prts(i)%collided)) then
!               prts(i)%x = 0.01D0*prts(i)%vel +prts(i)%x
!            else
!               prts(i)%collided = .false.
!            end if
!
!            print *, prts(i)%x,time
!         end do


      END MODULE PRT2MOD