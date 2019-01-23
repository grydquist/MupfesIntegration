!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     This routine contains a class for time integration. The algorithm
!     is based on Generalized-\alpha method that contains predictor,
!     initiator, and corrector routines. For usage and testing the class
!     see TEST_ITGMOD at the end of this file. 
!      
!--------------------------------------------------------------------
      MODULE ITGMOD
      USE VARMOD
      IMPLICIT NONE

!     For time integration of PDEs
      TYPE itgType
!        \alpha_f
         REAL(KIND=8) af
!        \alpha_m
         REAL(KIND=8) am
!        \beta
         REAL(KIND=8) beta
!        \gamma
         REAL(KIND=8) gam
!        \ro_{infinity}
         REAL(KIND=8) :: roInf = 0.2D0
      CONTAINS 
!        Predicting variables at next time step
         PROCEDURE :: predict => predictItg
!        Initiating current version of variables
         PROCEDURE :: initiate => initiateItg
!        Correcting variables at the next time 
         PROCEDURE :: correctSItg
         PROCEDURE :: correctVItg
         GENERIC :: correct => correctSItg, correctVItg
      END TYPE itgType

      INTERFACE itgType
         PROCEDURE :: newItg
      END INTERFACE itgType
!---------------------------------------------------------------------      
!     Number of time steps
      INTEGER nTS
!     Current time step
      INTEGER :: cTS = 0
!     Time step size
      REAL(KIND=8) dt
!     Time
      REAL(KIND=8) :: time = 0D0

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
      FUNCTION newItg(itgO, roInf) RESULT(itg)
      TYPE(itgType) :: itg
      INTEGER, INTENT(IN) :: itgO
      REAL(KIND=8), INTENT(IN), OPTIONAL :: roInf

      IF (PRESENT(roInf)) itg%roInf = roInf
      
      SELECT CASE (itgO)
      CASE(1)
!     For first order equations
         itg%am = 5D-1*(3D0 - itg%roInf)/(1D0 + itg%roInf)
      CASE(2)
!     For second order equations
         itg%am = (2D0 - itg%roInf)/(1D0 + itg%roInf)
      CASE DEFAULT
         io%e = "Integration scheme must be 1 or 2"
      END SELECT

      itg%af    = 1D0/(1D0 + itg%roInf)
      itg%beta  = 25D-2*(1D0 + itg%am - itg%af)**2D0
      itg%gam   = 5D-1 + itg%am - itg%af

      RETURN
      END FUNCTION newItg
!---------------------------------------------------------------------
!     This is the predictor      
      SUBROUTINE predictItg(itg,s)
      CLASS(itgType), INTENT(IN) :: itg
      CLASS(gVarType), INTENT(INOUT) :: s

      REAL(KIND=8) coef
        
      s%n%sP  = s%o%sP
      coef    = (itg%gam - 1D0)/itg%gam
      s%An%sP = s%Ao%sP*coef
      coef    = dt*dt*(5D-1*itg%gam - itg%beta)/(itg%gam - 1D0)
      s%Dn%sP = s%Do%sP + s%n%sP*dt + s%An%sP*coef

      RETURN
      END SUBROUTINE predictItg
!---------------------------------------------------------------------
!     This is the initiator      
      SUBROUTINE initiateItg(itg,s)
      CLASS(itgType), INTENT(IN) :: itg
      CLASS(gVarType), INTENT(INOUT) :: s
    
      s%A%sP = s%Ao%sP*(1D0 - itg%am) + s%An%sP*itg%am
      IF (LOWER(s%name) .EQ. 'pressure') THEN
         s%sP = s%n%sP
      ELSE
         s%sP = s%o%sP*(1D0 - itg%af) + s%n%sP*itg%af
      END IF
      s%D%sP = s%Do%sP*(1D0 - itg%af) + s%Dn%sP*itg%af

      RETURN
      END SUBROUTINE initiateItg
!---------------------------------------------------------------------
!     This is the corrector.
      SUBROUTINE correctSItg(itg, s, R)
      CLASS(itgType), INTENT(IN) :: itg
      CLASS(gVarType), INTENT(INOUT) :: s
      REAL(KIND=8), INTENT(IN) :: R(:)
      
      REAL(KIND=8) cY, cD

      IF (SIZE(R) .NE. SIZE(s%sP)) io%e = "correctSItg: Incompatible R"

      cY = itg%gam*dt
      cD = itg%beta*dt*dt
      
      s%An%sP = s%An%sP - R
      s%n%sP  = s%n%sP  - R*cY
      s%Dn%sP = s%Dn%sP - R*cD
  
      RETURN  
      END SUBROUTINE correctSItg
!---------------------------------------------------------------------
!     This is the corrector.
      SUBROUTINE correctVItg(itg, s, R)
      CLASS(itgType), INTENT(IN) :: itg
      CLASS(gVarType), INTENT(INOUT) :: s(:)
      REAL(KIND=8), INTENT(IN) :: R(:)
      
      INTEGER i, nV, nR, nNo, a, f, l
      REAL(KIND=8) cD
      INTEGER, ALLOCATABLE :: ds(:), fs(:), ls(:)
      REAL(KIND=8), ALLOCATABLE :: cY(:)

      nV = SIZE(s)
      IF (nV .LT. 1) RETURN

      ALLOCATE(cY(nV), ds(nV), fs(nV), ls(nV))
      cD  = itg%beta*dt*dt
      cY  = itg%gam*dt
      nNo = s(1)%dmn%dof
      a   = 0
      DO i=1, nV
         a     = a + SIZE(s(i)%sP)
         ds(i) = s(i)%dof
         IF (s(i)%dmn%dof.NE.nNo) io%e = "correctVItg: inconsistent dof"
         IF (LOWER(s(i)%name) .EQ. 'pressure') cY(i) = itg%af*cY(i)
      END DO
      nR = SIZE(R)
      IF (a .NE. nR) io%e = "correctVItg: Incompatible R"

      l  = 0
      ls = 0
      DO a=1, nNo
         fs = ls + 1
         ls = ls + ds
         DO i=1, nV
            f = l + 1
            l = l + ds(i)

            s(i)%An%sP(fs(i):ls(i)) = s(i)%An%sP(fs(i):ls(i)) 
     2         - R(f:l)
            s(i)% n%sP(fs(i):ls(i)) = s(i)% n%sP(fs(i):ls(i)) 
     2         - R(f:l)*cY(i)
            s(i)%Dn%sP(fs(i):ls(i)) = s(i)%Dn%sP(fs(i):ls(i)) 
     2         - R(f:l)*cD
         END DO
      END DO
 
      RETURN  
      END SUBROUTINE correctVItg
!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_ITGMOD(prs)
      TYPE(gVarType), INTENT(INOUT) :: prs

      TYPE(itgType) itg
      REAL(KIND=8), ALLOCATABLE :: R(:)

!     All the needed structures/classes
      itg = itgType(1) ! An integrator for first order PDE
      io%o = "newItg: "//CLR("(PASSED)",3)
      prs%Ao%sP = 0D0
      CALL itg%predict(prs)
      io%o = "itg%predictor: "//CLR("(PASSED)",3)
      CALL itg%initiate(prs)
      io%o = "itg%predictor: "//CLR("(PASSED)",3)
      ALLOCATE(R(prs%dmn%dof))
      R = 1D0
      CALL itg%correct(prs,R)
      IF (ANY(prs%An%sP.NE.-1D0)) io%e = "Issue with %corrector"
      io%o = "itg%predictor: "//CLR("(PASSED)",3)

      RETURN
      END SUBROUTINE TEST_ITGMOD

      END MODULE ITGMOD
