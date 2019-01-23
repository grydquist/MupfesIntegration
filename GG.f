!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     Module for prescribed boundary conditions. To be used for time
!     varying BC as a function of space and/or time. 
!     For usage and testing this class, see TEST_GGMOD at the end of 
!     this file. 
!      
!---------------------------------------------------------------------
      MODULE GGMOD
      USE VARMOD
      IMPLICIT NONE

!     Keeping track of time-dependent variables using Fourier transform
      TYPE gtType
         PRIVATE
!        Degrees of freedom
         INTEGER dof
!        Number of Fourier coefficient
         INTEGER :: n = 0
!        Initial value (offset)
         REAL(KIND=8), ALLOCATABLE :: qi(:)
!        Time derivative of linear part
         REAL(KIND=8), ALLOCATABLE :: qs(:)
!        Period
         REAL(KIND=8) T
!        Initial time
         REAL(KIND=8) ti
!        Imaginary part of coefficint
         REAL(KIND=8), ALLOCATABLE :: i(:,:)
!        Real part of coefficint
         REAL(KIND=8), ALLOCATABLE :: r(:,:)
      CONTAINS
!        Evaluates at a given time point
         PROCEDURE :: eval => evalGt
!        Deallocates all structures
         PROCEDURE :: free => freeGt
      END TYPE gtType

      INTERFACE gtType
         PROCEDURE :: newGt
      END INTERFACE gtType
!---------------------------------------------------------------------
!     For prescribed values that are a function of x and t. 
      TYPE gmType
         PRIVATE
!        Degrees of freedom of g    
         INTEGER dof
!        Number of time points to be read
         INTEGER :: nTP = 0
!        The period of data
         REAL(KIND=8) period
!        Time points
         REAL(KIND=8), ALLOCATABLE :: t(:)
!        Displacements at each direction, location, and time point
         REAL(KIND=8), ALLOCATABLE :: d(:,:,:)
      CONTAINS
!        Evaluates at a given time point
         PROCEDURE :: eval => evalGm
!        Deallocates all structures
         PROCEDURE :: free => freeGm
      END TYPE gmType

      INTERFACE gmType
         PROCEDURE :: newGm
      END INTERFACE gmType
!---------------------------------------------------------------------
!     For prescribed values that are a function of x and t. 
      TYPE gxType
         PRIVATE
!        Degrees of freedom of g
         INTEGER :: dof
!        Values of at nodal points 
         REAL(KIND=8), ALLOCATABLE :: s(:,:)
      CONTAINS
!        Deallocates all structures
         PROCEDURE :: free => freeGx
      END TYPE gxType

      INTERFACE gxType
         PROCEDURE :: newGx
      END INTERFACE gxType
!---------------------------------------------------------------------
!     General prescribed BC that acts as a wrapper for gx, gt, and gm
      TYPE ggType
         PRIVATE
!        Face associated with this BC
         LOGICAL :: mapped = .FALSE.
!        Face associated with this BC
         INTEGER iFa
!        Mesh associated with this BC
         INTEGER iM
!        Constant (steady) value
         REAL(KIND=8), ALLOCATABLE, PUBLIC :: gs(:)
!        Spatial dependent: gx=gx(x)
         TYPE(gxType), ALLOCATABLE :: gx
!        Time dependent: gt=gt(t)
         TYPE(gtType), ALLOCATABLE :: gt
!        Spatio-temporal dependent: gm=gm(x,t)
         TYPE(gmType), ALLOCATABLE :: gm
!        Pointer to the domain associated with this BC
         TYPE(dmnType), POINTER :: dmn => NULL()
      CONTAINS 
!        Computes the profile
         PROCEDURE :: cProfile => cProfileGg
!        Provides the prescribed value at x and t
         PROCEDURE :: eval => evalGg
!        Deallocates all structures
         PROCEDURE :: free => freeGg
!        Maps variables from nNo to dof
         PROCEDURE, PRIVATE :: mapSGg
         PROCEDURE, PRIVATE :: mapVGg
         PROCEDURE :: map => mapGg
      END TYPE ggType

      INTERFACE ggType
         PROCEDURE :: newGg, newGgFl
      END INTERFACE ggType

      CONTAINS 
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newGt(fid) RESULT(gt)
      INTEGER, INTENT(IN) :: fid
      TYPE(gtType) :: gt

      INTEGER i, n, np
      REAL(KIND=8) tmp, kn, ko
      REAL(KIND=8), ALLOCATABLE :: t(:), q(:,:), s(:)
      
      READ(fid,*) i, np
      IF (np .LT. 2) io%e = "Enter dof,nPnts followed by nPts*(t Q)"
      IF (i.LT.1 .OR. i.GT.nsd) io%e = "newGt: Degrees of freedom"//
     2   " must be between 0 and "//nsd

      gt%n   = np
      gt%dof = i
      ALLOCATE(gt%r(i,np), gt%i(i,np), gt%qi(i), gt%qs(i), t(np), 
     2   q(i,np), s(i))

      DO i=1, np
         READ (fid,*) t(i), q(:,i)
      END DO
      CLOSE(fid)

      gt%ti = t(1)
      gt%T  = t(np) - t(1)
      gt%qi = q(:,1)
      gt%qs = (q(:,np) - q(:,1))/gt%T

      DO i=1, np
         t(i)   = t(i) - gt%ti
         q(:,i) = q(:,i) - gt%qi - gt%qs*t(i)
      END DO

      DO n=1, gt%n
         tmp = REAL(n-1,8)
         gt%r(:,n) = 0D0
         gt%i(:,n) = 0D0
         DO i=1, np-1
            ko = 2D0*pi*tmp*t(i)/gt%T
            kn = 2D0*pi*tmp*t(i+1)/gt%T
            s  = (q(:,i+1) - q(:,i))/(t(i+1) - t(i))

            IF (n .EQ. 1) THEN
               gt%r(:,n) = gt%r(:,n) + 5D-1*(t(i+1) - t(i))*
     2            (q(:,i+1) + q(:,i))
            ELSE
               gt%r(:,n) = gt%r(:,n) + s*(COS(kn) - COS(ko))
               gt%i(:,n) = gt%i(:,n) - s*(SIN(kn) - SIN(ko))
            END IF
         END DO

         IF (n .EQ. 1) THEN
            gt%r(:,n) = gt%r(:,n)/gt%T
         ELSE
            gt%r(:,n) = 5D-1*gt%r(:,n)*gt%T/(pi*pi*tmp*tmp)
            gt%i(:,n) = 5D-1*gt%i(:,n)*gt%T/(pi*pi*tmp*tmp)
         END IF
      END DO
      DEALLOCATE(s, q, t)

      RETURN
      END FUNCTION newGt
!---------------------------------------------------------------------
!     This is to calculate flow rate and flow acceleration (IFFT)
      PURE SUBROUTINE evalGt(gt, time, Y, A)
      CLASS(gtType), INTENT(IN) :: gt
      REAL(KIND=8), INTENT(IN) :: time
      REAL(KIND=8), INTENT(OUT) :: Y(gt%dof), A(gt%dof)

      INTEGER i
      REAL(KIND=8) t, tmp, K, kd

      t   = DMOD(time - gt%ti, gt%T)
      tmp = 2D0*pi/gt%T
      Y   = gt%qi + t*gt%qs
      A   = gt%qs
      DO i=1, gt%n
         kd = tmp*REAL(i-1,8)
         K  = t*kd
         Y  = Y +  gt%r(:,i)*COS(K) - gt%i(:,i)*SIN(K)
         A  = A - (gt%r(:,i)*SIN(K) + gt%i(:,i)*COS(K))*kd
      END DO

      RETURN 
      END SUBROUTINE evalGt
!---------------------------------------------------------------------
!     Communicating time-depandant BC data
      SUBROUTINE freeGt(gt)
      CLASS(gtType), INTENT(INOUT) :: gt
      
      IF (ALLOCATED(gt%r))  DEALLOCATE(gt%r)
      IF (ALLOCATED(gt%i))  DEALLOCATE(gt%i)
      IF (ALLOCATED(gt%qi)) DEALLOCATE(gt%qi)
      IF (ALLOCATED(gt%qs)) DEALLOCATE(gt%qs)
      gt%n = 0

      RETURN
      END SUBROUTINE freeGt

!#####################################################################
!     To read mb structure from a file
      FUNCTION newGm(lM, iFa, fid) RESULT(gm)
      TYPE(gmType) :: gm
      INTEGER, INTENT(IN) :: fid
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: iFa

      INTEGER i, a, Ac, nNo
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp
      INTEGER, ALLOCATABLE :: gptr(:), lptr(:)
      REAL(KIND=8), ALLOCATABLE :: tmp(:,:,:)
      ASSOCIATE(fa => lM%fa(iFa))

      READ (fid,*) gm%dof, gm%nTP, nNo
      ALLOCATE(gm%t(gm%nTP))
!     I am seting all the nodes to zero just in case a node is not set
      DO i=1, gm%nTP
         READ (fid,*) rtmp
         gm%t(i) = rtmp
         IF (i .EQ. 1) THEN
            IF (.NOT.ISZERO(rtmp)) io%e = "First time step"//
     2         " should be zero in <"//TRIM(ctmp)//">"
         ELSE
            rtmp = rtmp - gm%t(i-1)
            IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) io%e = "Non-in"//
     2         "creasing time trend is found in <"//TRIM(ctmp)//">"
         END IF
      END DO
      gm%period = gm%t(gm%nTP)

!     Preparing the pointer array
      IF (cm%mas()) THEN
         IF (gm%dof.LT.1 .OR. gm%dof.GT.nsd) io%e = "newGm: "//
     2      "Degrees of freedom must be between 0 and "//nsd
         ALLOCATE(gptr(lM%gnNo))         
         gptr = 0
         DO a=1, nNo
            READ(fid,*) Ac
            IF (Ac.GT.lM%gnNo .OR. Ac.LE.0) io%e = "Entry "//a//
     2         " is out of bound in "//ctmp
            gptr(Ac) = a
            DO i=1, gm%nTP
               READ(fid,*) 
            END DO
         END DO
!     Going back to the begining of the block
         REWIND(fid)
         DO i=1, gm%nTP + 1
            READ(fid,*)
         END DO
      ELSE
         ALLOCATE(gptr(0))
      END IF
      ALLOCATE(lptr(lM%nNo))
      CALL cm%local(lptr, gptr, lM%gN)
!     Reading all the data into tmp
      ALLOCATE(gm%d(gm%dof,fa%nNo,gm%nTP), tmp(gm%dof,nNo,gm%nTP))
      tmp = 0D0
      DO a=1, nNo
         READ(fid,*) Ac
         DO i=1, gm%nTP
            READ (fid,*) tmp(:,a,i)
         END DO
      END DO
      CLOSE(fid)
!     Assigning corresponding entries to %d
      DO a=1, fa%nNo
         Ac = lptr(fa%gN(a))
         gm%d(:,a,:) = tmp(:,Ac,:)
      END DO
      DEALLOCATE(gptr, lptr, tmp)

      RETURN
      END ASSOCIATE
      END FUNCTION newGm
!---------------------------------------------------------------------
!     Calculating boundary values from gm data structure
      PURE SUBROUTINE evalGm(gm, time, Y, A)
      CLASS(gmType), INTENT(IN) :: gm
      REAL(KIND=8), INTENT(IN) :: time
      REAL(KIND=8), INTENT(OUT) :: Y(gm%dof,SIZE(gm%d,2)),
     2   A(gm%dof,SIZE(gm%d,2))

      INTEGER b, i, nNo
      REAL(KIND=8) t, tmp, delT

      t = DMOD(time,gm%period)
      DO i=1, gm%nTP - 1
         IF (gm%t(i+1) .GE. t) THEN
            Y = 0D0
            A = 0D0
            EXIT
         END IF
      END DO
      delT = gm%t(i+1) - gm%t(i)
      tmp  = (t - gm%t(i))/delT
      nNo  = SIZE(gm%d,2)
      DO b=1, nNo
         Y(:,b) = tmp*gm%d(:,b,i+1) + gm%d(:,b,i)*(1D0-tmp)
         A(:,b) =    (gm%d(:,b,i+1) - gm%d(:,b,i))/delT
      END DO

      RETURN
      END SUBROUTINE evalGm
!---------------------------------------------------------------------
!     Communicating time-depandant BC data
      SUBROUTINE freeGm(gm)
      CLASS(gmType), INTENT(INOUT) :: gm
      
      IF (ALLOCATED(gm%t)) DEALLOCATE(gm%t)
      IF (ALLOCATED(gm%d)) DEALLOCATE(gm%d)
      gm%nTP = 0

      RETURN
      END SUBROUTINE freeGm

!#####################################################################

      FUNCTION newGx(lM, iFa, fid) RESULT(gx)
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: iFa, fid
      TYPE(gxType) :: gx
 
      INTEGER a, Ac, nNo
      REAL(KIND=8), ALLOCATABLE :: ltmp(:,:), tmp(:), gtmp(:,:)
      ASSOCIATE(fa => lM%fa(iFa))

      READ(fid,*) gx%dof, nNo
      IF (cm%mas()) THEN
         IF (gx%dof.LT.1 .OR. gx%dof.GT.nsd) io%e = "newGx: "//
     2      "Degrees of freedom must be between 0 and "//nsd
         ALLOCATE(gtmp(gx%dof,lM%gnNo), tmp(gx%dof))
         gtmp = 0D0
         DO a=1, nNo
            READ(fid,*) Ac, tmp
            IF (Ac.GT.lM%gnNo .OR. Ac.LE.0) io%e = "newGx: entry "//a
     2         //" is out of bound"
            gtmp(:,Ac) = tmp
         END DO
         DEALLOCATE(tmp)
      ELSE
         ALLOCATE(gtmp(0,0))
      END IF
      CLOSE(fid)
      ALLOCATE(ltmp(gx%dof,lM%nNo))
      CALL cm%local(ltmp, gtmp, lM%gN)

      ALLOCATE(gx%s(gx%dof,fa%nNo))
      DO a=1, fa%nNo
         Ac = fa%gN(a)
         gx%s(:,a) = ltmp(:,Ac)
      END DO
      DEALLOCATE(ltmp, gtmp)

      RETURN
      END ASSOCIATE
      END FUNCTION newGx
!---------------------------------------------------------------------
      SUBROUTINE freeGx(gx)
      CLASS(gxType), INTENT(INOUT) :: gx
      
      IF (ALLOCATED(gx%s)) DEALLOCATE(gx%s)

      RETURN
      END SUBROUTINE freeGx

!#####################################################################

      FUNCTION newGg(dmn, fName, gs, fGx, fGt, fGm) RESULT(gg)
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      REAL(KIND=8), INTENT(IN) :: gs(:)
      CHARACTER(LEN=*), INTENT(IN) :: fName
      INTEGER, INTENT(IN) :: fGx, fGt, fGm
      TYPE(ggType) :: gg

      LOGICAL flag

      gg%dmn => dmn
      CALL dmn%find(fName, gg%iM, gg%iFa)
      ALLOCATE(gg%gs, SOURCE=gs)

      INQUIRE(UNIT=fGm,OPENED=flag)
      IF (flag) THEN
         ALLOCATE(gg%gm)
         gg%gm = gmType(dmn%msh(gg%iM), gg%iFa, fGm)
      END IF

      INQUIRE(UNIT=fGt,OPENED=flag)
      IF (flag) THEN
         ALLOCATE(gg%gt)
         gg%gt = gtType(fGt)
      END IF

      INQUIRE(UNIT=fGx,OPENED=flag)
      IF (flag) THEN
         ALLOCATE(gg%gx)
         gg%gx = gxType(dmn%msh(gg%iM), gg%iFa, fGx)
      END IF

      RETURN
      END FUNCTION newGg
!---------------------------------------------------------------------
      FUNCTION newGgFl(lst, dmn, fName) RESULT(gg)
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      CHARACTER(LEN=*), INTENT(IN) :: fName
      TYPE(ggType) :: gg
      
      LOGICAL zp, flx
      INTEGER i
      CHARACTER(LEN=stdL) ctmp, stmp
      TYPE(lstType), POINTER :: lPtr
      TYPE(fileType) f

      gg%dmn => dmn
      CALL dmn%find(fName,gg%iM,gg%iFa)
      ctmp   = 'steady'
      lPtr  => lst%get(ctmp,"Time dependency")
      SELECT CASE (ctmp)
      CASE ('steady')
         lPtr => lst%get(f,"Value",1)
         stmp = f%name()
         i    = checkNoNumbers(stmp)
         ALLOCATE(gg%gs(i)) 
         READ(stmp,*) gg%gs
      CASE ('coupled')
         ALLOCATE(gg%gs(1)) 
      CASE ('unsteady')
         lPtr => lst%get(f,"Temporal values file path",1)
         CALL f%open('r')
         ALLOCATE(gg%gt)
         gg%gt = gtType(f%id())
      CASE ('general')
         lPtr =>lst%get(f,"Temporal and spatial values file path",1)
         CALL f%open('r')
         ALLOCATE(gg%gm)
         gg%gm = gmType(dmn%msh(gg%iM), gg%iFa, f%id())
      CASE DEFAULT
         io%e=TRIM(lst%ping("Time dependency",lPtr))//" Unexpected type"
      END SELECT

!     To zero-out perimeter or not. Default is .true. for Dir 
      zp = .FALSE.
      lPtr => lst%get(ctmp,"type")
      IF (ctmp.EQ."dirichlet" .OR. ctmp.EQ."dir") zp = .TRUE.
      lPtr => lst%get(zp,"Zero out perimeter")

!     To impose value or flux
      flx = .FALSE.
      lPtr => lst%get(flx,"Impose flux")

!     Reading the spatial profile: flat/para/ud
      ctmp = "flat"
      lPtr => lst%get(ctmp,"Profile")
      SELECT CASE (ctmp)
      CASE ('flat')
         CALL gg%cProfile(.FALSE.,zp,flx)
      CASE ('parabolic') 
         CALL gg%cProfile(.TRUE.,zp,flx)
      CASE ('user_defined') 
         ALLOCATE(gg%gx)
         lPtr => lst%get(f,"Spatial profile file path",1)
         CALL f%open('r')
         gg%gx = gxType(dmn%msh(gg%iM), gg%iFa, f%id())
      CASE DEFAULT
         io%e = TRIM(lst%ping("Profile",lPtr))//" Unexpected profile"
      END SELECT

      RETURN
      END FUNCTION newGgFl
!---------------------------------------------------------------------
!     This is to calculate flow rate and flow acceleration
      SUBROUTINE evalGg(gg, time, Y, A)
      CLASS(ggType), INTENT(IN) :: gg
      REAL(KIND=8), INTENT(IN) :: time
      REAL(KIND=8), INTENT(INOUT) :: Y(:,:)
      REAL(KIND=8), INTENT(INOUT), OPTIONAL, TARGET :: A(:,:)

      LOGICAL flag
      INTEGER b, Ac, dof, ldof, tdof, xdof, lnNo
      REAL(KIND=8), ALLOCATABLE :: Yl(:,:), Al(:,:), Yt(:), At(:),
     2   Yv(:,:), Av(:,:)
      ASSOCIATE(lM => gg%dmn%msh(gg%iM), fa => lM%fa(gg%iFa))

      flag = PRESENT(A)
      IF (flag) THEN
         IF (SIZE(A,1).NE.SIZE(Y,1) .OR. SIZE(A,2).NE.SIZE(Y,2)) 
     2      io%e = "evalGg: Unexpected SIZE(A)"
      END IF
      dof = SIZE(Y,1)

!     DOF of Yl/Al should be either equal to dof of Y/A computed above
!     or should be 1, which in that case a normal vector contracted with
!     Yl/Al to bring its dof up to Y/A. 

      lnNo = fa%nNo
      IF (gg%mapped) lnNo = fa%dof
      IF (ALLOCATED(gg%gm)) THEN ! The general case
         ldof = gg%gm%dof
         ALLOCATE(Yl(ldof,lnNo), Al(ldof,lnNo))
         CALL gg%gm%eval(time, Yl, Al)
         RETURN
      ELSE
         IF (ALLOCATED(gg%gt)) THEN ! TIme dependency
            tdof = gg%gt%dof
            ALLOCATE(Yt(tdof), At(tdof))
            CALL gg%gt%eval(time, Yt, At)
         ELSE ! steady
            tdof = SIZE(gg%gs)
            ALLOCATE(Yt(tdof), At(tdof))
            Yt = gg%gs
            At = 0D0
         END IF
         xdof = 1
         IF (ALLOCATED(gg%gx)) xdof = gg%gx%dof
!     xdof should be compatible with tdof in the following sense
         IF (xdof.NE.1 .AND. tdof.NE.1 .AND. xdof.NE.tdof) io%e =
     2      "evalGG: Time and Spatial profile have incompatible DOF"

         ldof = MAX(xdof,tdof)
         ALLOCATE(Yl(ldof,lnNo), Al(ldof,lnNo))
         IF (ALLOCATED(gg%gx)) THEN ! Spatial dependency
            DO b=1, lnNo
               Yl(:,b) = Yt*gg%gx%s(:,b)
               Al(:,b) = At*gg%gx%s(:,b)
            END DO
         ELSE 
            DO b=1, lnNo
               Yl(:,b) = Yt
               Al(:,b) = At
            END DO
         END IF
         DEALLOCATE(Yt, At)
      END IF

!     dof and ldof should be compatible in the following sense
      IF (ldof.EQ.1 .AND. dof.EQ.nsd) THEN
         ALLOCATE(Yv(nsd,lnNo), Av(nsd,lnNo))
         DO b=1, lnNo
            Yv(:,b) = Yl(1,b)*fa%nV(:,b)
            Av(:,b) = Al(1,b)*fa%nV(:,b)
         END DO
         DEALLOCATE(Yl, Al)
      ELSE IF (ldof.EQ.dof) THEN
         CALL MOVE_ALLOC(Yl,Yv)
         CALL MOVE_ALLOC(Al,Av)
      ELSE
         io%e = "evalGG: variable and imposed BC have incompatible DOF"
      END IF

      IF (gg%mapped) THEN
         IF (SIZE(Y,2) .NE. fa%dof) io%e = "ecalGg: ~fa%dof for mapped"
         Y = Yv
         IF (flag) A = Av
      ELSE
         IF (SIZE(Y,2) .NE. lM%nNo) io%e = "ecalGg: ~lM%nNo for n.mappd"
         DO b=1, fa%nNo
            Ac = fa%gN(b)
            Y(:,Ac) = Yv(:,b)
            IF (flag) A(:,Ac) = Av(:,b)
         END DO
      END IF
      DEALLOCATE(Yv, Av)

      RETURN 
      END ASSOCIATE
      END SUBROUTINE evalGg
!---------------------------------------------------------------------
!     Compute a profile provided whether the profile should be flat or
!     parabolic, the primeter should be zero-out or not, and the flux 
!     should be imposed or the absolute value.
      SUBROUTINE cProfileGg(gg, parabolic, zeroP, iFlx)
      CLASS(ggType), INTENT(INOUT) :: gg
      LOGICAL, INTENT(IN), OPTIONAL :: parabolic, zeroP, iFlx
      
      LOGICAL fParabolic, fZeroP, fIFlx
      INTEGER iFa, i, a, b, Ac, j, nNo
      REAL(KIND=8) tmp, nV(nsd), center(nsd), maxN
      LOGICAL, ALLOCATABLE :: gN(:)
      REAL(KIND=8), ALLOCATABLE :: s(:), sV(:,:), sVl(:,:)
      ASSOCIATE(lM => gg%dmn%msh(gg%iM), fa => lM%fa(gg%iFa))

      nNo = lM%nNo
      IF (ALLOCATED(gg%gx)) CALL gg%gx%free()
      IF (.NOT.ALLOCATED(gg%gx)) ALLOCATE(gg%gx)
      gg%gx%dof = 1
      ALLOCATE(gg%gx%s(gg%gx%dof,fa%nNo))

      fParabolic = .TRUE.
      fZeroP     = .TRUE.
      fIFlx      = .FALSE.
      IF (PRESENT(parabolic)) fParabolic = parabolic
      IF (PRESENT(zeroP))     fZeroP     = zeroP
      IF (PRESENT(iFlx))      fIFlx      = iFlx

      ALLOCATE(s(nNo))
      s = 0D0
!     Flat profile
      IF (.NOT.fParabolic) THEN
         DO a=1, fa%nNo
            Ac    = fa%gN(a)
            s(Ac) = 1D0
         END DO 
      ELSE
!     Here is the method that is used for imposing parabolic profile:
!     1- Find the coordinate of the points on the boundary 2- find unit 
!     vector from center to each of points, i, on the boundary: ew(i)
!     3- maximize ew(i).e where e is the unit vector from current 
!     point to the center 4- Use the point i as the diam here
         DO i=1, nsd
            center(i) = lM%integ(gg%iFa, lM%x, i)/fa%area
         END DO
         center = cm%reduce(center)
         ALLOCATE(gN(nNo), sVl(nsd,fa%nNo))
!     gN is true if a node located on the boundary (beside fa)
         gN = .FALSE.
         gN(fa%ring) = .TRUE.

!     "j" is a counter for the number of nodes that are located on the 
!     boundary of fa and sVl contains the list of their coordinates
         j   = 0
         sVl = 0D0
         DO a=1, fa%nNo
            Ac = fa%gN(a)
            IF (gN(Ac)) THEN
               j        = j + 1
               sVl(:,j) = lM%x(:,Ac)
            END IF
         END DO
         IF (.NOT.cm%seq()) THEN
            sV = cm%gather(sVl(:,1:j),toAll=.TRUE.)
            j  = SIZE(sV(1,:))
            DEALLOCATE(sVl)
         ELSE
            CALL MOVE_ALLOC(sVl,sV)
         END IF
         IF (j .EQ. 0) io%e = "No perimeter found for face: "//fa%name

!     Computing the vector from center to perimeter
         DO a=1, j
            sV(:,a) = sV(:,a) - center
         END DO
!     "s" is going to keep the ew.e value
         DO a=1, fa%nNo
            Ac = fa%gN(a)
            nV = lM%x(:,Ac) - center
            maxN = -HUGE(maxN)
            DO b=1, j
               tmp = DOT_PRODUCT(nV,sV(:,b))/NORM2(sV(:,b))
               IF (tmp .GT. maxN) THEN
                  maxN = tmp
                  i = b
               END IF
            END DO
            s(Ac) = MAX(0D0,1D0 - (NORM2(nV)/NORM2(sV(:,i)))**2D0)
         END DO
         DEALLOCATE(sV)
      END IF

!     Now correcting the inlet BC for the inlet ring
      IF (fZeroP) s(fa%ring) = 0D0

!     Normalizing the profile for flux
      tmp = 1D0
      IF (fIFlx) THEN
         tmp = lM%integ(gg%iFa,s)
         tmp = cm%reduce(tmp)
         IF (ISZERO(tmp)) THEN
            tmp = 1D0
            io%w = "Using face <"//TRIM(fa%name)//
     2         "> to impose BC led to no non-zero node."
         END IF
      END IF
      DO a=1, fa%nNo
         Ac = fa%gN(a)
         gg%gx%s(:,a) = s(Ac)/tmp
      END DO

      RETURN
      END ASSOCIATE
      END SUBROUTINE cProfileGg
!---------------------------------------------------------------------
!     Communicating time-depandant BC data
      SUBROUTINE freeGg(gg)
      CLASS(ggType), INTENT(INOUT) :: gg

      IF (ALLOCATED(gg%gm)) THEN
         CALL gg%gm%free()
         DEALLOCATE(gg%gm)
      END IF
      IF (ALLOCATED(gg%gt)) THEN
         CALL gg%gt%free()
         DEALLOCATE(gg%gt)
      END IF
      IF (ALLOCATED(gg%gx)) THEN
         CALL gg%gx%free()
         DEALLOCATE(gg%gx)
      END IF

      IF (ALLOCATED(gg%gs)) DEALLOCATE(gg%gs)
      
      gg%dmn => NULL()
      gg%mapped = .FALSE.

      RETURN
      END SUBROUTINE freeGg
!---------------------------------------------------------------------
!     Distributing data to neighbouring partitions
      SUBROUTINE mapSGg(gg, s)
      CLASS(ggType), INTENT(INOUT) :: gg
      REAL(KIND=8), INTENT(INOUT), ALLOCATABLE :: s(:,:)

      INTEGER a, dof
      REAL(KIND=8), ALLOCATABLE :: sG(:,:)
      INTEGER, ALLOCATABLE :: map(:)
      ASSOCIATE(lM => gg%dmn%msh(gg%iM), fa => lM%fa(gg%iFa))

      dof = SIZE(s,1)
      ALLOCATE(map(fa%nNo))
      IF (cm%mas()) THEN
         ALLOCATE(sG(dof,gg%dmn%gdof))
         sG = 0D0
      ELSE
         ALLOCATE(sG(0,0))
      END IF
      DO a=1, fa%nNo
         map(a) = gg%dmn%guP(fa%uP(a))
      END DO
      CALL cm%global(s, sG, map)
      DEALLOCATE(s, map)
      ALLOCATE(s(dof,fa%dof), map(fa%dof))
      DO a=1, fa%dof
         map(a) = gg%dmn%guP(fa%uP(a))
      END DO
      CALL cm%local(s, sG, map)
      DEALLOCATE(sG, map)

      RETURN
      END ASSOCIATE
      END SUBROUTINE mapSGg
!---------------------------------------------------------------------
!     Distributing data to neighbouring partitions
      SUBROUTINE mapVGg(gg, s)
      CLASS(ggType), INTENT(INOUT) :: gg
      REAL(KIND=8), INTENT(INOUT), ALLOCATABLE :: s(:,:,:)

      INTEGER i, m, n, dof, nNo
      REAL(KIND=8), ALLOCATABLE :: tmp(:,:), stmp(:,:,:)

      m   = SIZE(s,1)
      n   = SIZE(s,3)
      nNo = gg%dmn%msh(gg%iM)%fa(gg%iFa)%nNo
      dof = gg%dmn%msh(gg%iM)%fa(gg%iFa)%dof
      CALL MOVE_ALLOC(s,stmp)
      ALLOCATE(s(m,dof,n))
      DO i=1, n
         ALLOCATE(tmp(m,nNo))
         tmp = stmp(:,:,i)
         CALL mapSGg(gg, tmp)
         s(:,:,i) = tmp
         DEALLOCATE(tmp)
      END DO
      DEALLOCATE(stmp)

      RETURN
      END SUBROUTINE mapVGg
!---------------------------------------------------------------------
!     Distributing data to neighbouring partitions
      SUBROUTINE mapGg(gg)
      CLASS(ggType), INTENT(INOUT) :: gg

      IF (gg%mapped) io%e = "mapGg: already mapped"
      IF (ALLOCATED(gg%gm)) CALL mapVGg(gg, gg%gm%d)
      IF (ALLOCATED(gg%gx)) CALL mapSGg(gg, gg%gx%s)
      gg%mapped = .TRUE.

      RETURN
      END SUBROUTINE mapGg

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_GGMOD(dmn, gg)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(ggType), INTENT(OUT) :: gg

      INTEGER nNo
      REAL(KIND=8) gs(1)
      REAL(KIND=8), ALLOCATABLE :: Y(:,:), A(:,:)

      gs = (/1D-1/)
      gg = ggType(dmn,'Face_1',gs,0,0,0) !., ., gs, fGx, fGt, fGm
      CALL gg%map() ! This is Dir BC
      IF (ANY(gg%gs .NE. gs)) io%e = "Issue with newGg"
      io%o = "newGg: "//CLR("(PASSED)",3)
      nNo = dmn%msh(gg%iM)%fa(gg%iFa)%nNo
      ALLOCATE(Y(1,nNo),A(1,nNo))
      CALL gg%eval(1D-1,Y,A) ! time, Y, A
      IF (ANY(Y(:,1).NE.gs) .OR. ANY(A.NE.0D0)) io%e = "Issue with eval"
      io%o = "gg%eval: "//CLR("(PASSED)",3)
      CALL gg%cProfile(.FALSE.,.FALSE.,.TRUE.) ! Flat, no-zp, impose_f
      IF (.NOT.ISZERO(gg%gx%s(1,:),1D0)) io%e = "Issue with cProfile"
      io%o = "gg%cProfile: "//CLR("(PASSED)",3)
     
      RETURN
      END SUBROUTINE TEST_GGMOD

      END MODULE GGMOD
