!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     eq program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     For stoing accessing material properties of different substances,
!     i.e. fluid, solid, and particle.
!     For usage and testing this class see TEST_MATMOD at the end of
!     this file. 
!      
!---------------------------------------------------------------------
      MODULE MATMOD
      USE LSTMOD
      IMPLICIT NONE

!     The total number of properties
      INTEGER, PARAMETER, PRIVATE :: mat_nPrp = 14

!     The identifier for various properties
      INTEGER, PARAMETER, PRIVATE :: mat_rho = 1, mat_mu = 2, 
     2   mat_A = 3, mat_fx = 4, mat_fy = 5, mat_fz = 6, 
     3   mat_K = 7, mat_Cp = 8, mat_q = 9, mat_beta = 10, 
     4   mat_E = 11, mat_nu = 12, mat_eta = 13, mat_D = 14

!     Whether a property has lower bound of 0D0
      LOGICAL, PARAMETER, PRIVATE :: mat_lb(mat_nPrp) = (/
     2   .TRUE., .TRUE., .FALSE., .FALSE., .FALSE., .FALSE., .TRUE., 
     3   .TRUE., .FALSE., .FALSE., .TRUE., .TRUE., .FALSE., .TRUE./)

!     The name associated with each properties
      CHARACTER(LEN=*), PARAMETER, PRIVATE :: mat_name(mat_nPrp) = (/
     2   "density",
     3   "dynamic viscosity",
     4   "ns source",
     5   "body force_x",
     6   "body force_y",
     7   "body force_z",
     8   "conductivity",
     9   "specific heat",
     z   "heat source",
     1   "backflow coefficient",
     2   "elasticity modulus",
     3   "poisson's ratio",
     4   "damping coefficient",
     6   "diameter"/)

!     Abstract material class
      TYPE matType
         PRIVATE
         CHARACTER(LEN=stdL) :: name = DEFAULT_NAME
         REAL(KIND=8) :: prp(mat_nPrp) = 0D0
      CONTAINS 
!        Density
         PROCEDURE :: rho => rhoMat
!        Dynamic viscosity
         PROCEDURE :: mu => muMat
!        Source term
         PROCEDURE :: A  => AMat
!        Body force_x
         PROCEDURE :: fx => fxMat
!        Body force_x
         PROCEDURE :: fy => fyMat
!        Body force_x
         PROCEDURE :: fz => fzMat
!        Conductivity
         PROCEDURE :: K => KMat
!        Specific heat
         PROCEDURE :: Cp => CpMat
!        Heat source term
         PROCEDURE :: q => qMat
!        Backflow stabilization
         PROCEDURE :: beta => betaMat
!        Elasticity modulus
         PROCEDURE :: E => EMat
!        Poisson's ratio
         PROCEDURE :: nu => nuMat
!        Damping coefficient
         PROCEDURE :: eta => etaMat
!        Diameter
         PROCEDURE :: D => DMat
!        Sets the value of a parameter provided its value and keyword
         PROCEDURE :: set => setMat
      END TYPE matType

      INTERFACE matType
         PROCEDURE :: newMat, newMatFl
      END INTERFACE matType
!---------------------------------------------------------------------
!     This is the material library
      TYPE(matType), ALLOCATABLE, TARGET, PRIVATE :: matLib(:)
!     FIND_MAT: To finds a material in this library provided its name
!     ADD_MAT: Adds a material to the library

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newMat(mName, pNames, pVals) RESULT(mat)
      CHARACTER(LEN=*), INTENT(IN) :: mName
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: pNames(:)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: pVals(:)
      TYPE(matType) :: mat
 
      INTEGER iPrp

      mat%name = LOWER(mName)
      IF (.NOT.PRESENT(pNames)) RETURN

      IF (.NOT.PRESENT(pVals)) io%e = "newMat: Both pNames and"//
     2   " pVals must be present at the same time"
      IF (SIZE(pNames) .NE. SIZE(pVals)) io%e = "newMat: "//
     2   "pNames and pVals must be of the same size"

      DO iPrp=1, SIZE(pNames)
         CALL mat%set(pNames(iPrp), pVals(iPrp))
      END DO

      RETURN
      END FUNCTION newMat
!---------------------------------------------------------------------
      FUNCTION newMatFl(lst, mName) RESULT(mat)
      TYPE(lstType), INTENT(INOUT) :: lst
      CHARACTER(LEN=*), INTENT(IN) :: mName
      TYPE(matType) :: mat
 
      INTEGER iPrp
      TYPE(lstType), POINTER :: lPtr

      mat%name = mName
      DO iPrp=1, mat_nPrp
         IF (mat_lb(iPrp)) THEN
            lPtr => lst%get(mat%prp(iPrp),mat_name(iPrp),lb=0D0)
         ELSE
            lPtr => lst%get(mat%prp(iPrp),mat_name(iPrp))
         END IF
      END DO

      RETURN
      END FUNCTION newMatFl
!---------------------------------------------------------------------
      PURE FUNCTION rhoMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_rho)
      END FUNCTION rhoMat
!---------------------------------------------------------------------
      PURE FUNCTION muMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_mu)
      END FUNCTION muMat
!---------------------------------------------------------------------
      PURE FUNCTION AMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_A)
      END FUNCTION AMat
!---------------------------------------------------------------------
      PURE FUNCTION fxMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_fx)
      END FUNCTION fxMat
!---------------------------------------------------------------------
      PURE FUNCTION fyMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_fy)
      END FUNCTION fyMat
!---------------------------------------------------------------------
      PURE FUNCTION fzMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_fz)
      END FUNCTION fzMat
!---------------------------------------------------------------------
      PURE FUNCTION KMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_K)
      END FUNCTION KMat
!---------------------------------------------------------------------
      PURE FUNCTION CpMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_Cp)
      END FUNCTION CpMat
!---------------------------------------------------------------------
      PURE FUNCTION qMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_q)
      END FUNCTION qMat
!---------------------------------------------------------------------
      PURE FUNCTION betaMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_beta)
      END FUNCTION betaMat
!---------------------------------------------------------------------
      PURE FUNCTION EMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_E)
      END FUNCTION EMat
!---------------------------------------------------------------------
      PURE FUNCTION nuMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_nu)
      END FUNCTION nuMat
!---------------------------------------------------------------------
      PURE FUNCTION etaMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_eta)
      END FUNCTION etaMat
!---------------------------------------------------------------------
      PURE FUNCTION DMat(mat) RESULT(s)
      CLASS(matType), INTENT(IN) :: mat
      REAL(KIND=8) s
      s = mat%prp(mat_D)
      END FUNCTION DMat
!---------------------------------------------------------------------
!     Set properties identified by kwd to u. If successful, TRUE is
!     returned, otherwise FALSE is returned
      SUBROUTINE setMat(mat, kwd, u)
      CLASS(matType), INTENT(INOUT) :: mat
      CHARACTER(LEN=*), INTENT(IN) :: kwd
      REAL(KIND=8), INTENT(IN) :: u
 
      INTEGER iPrp

      DO iPrp=1, mat_nPrp
         IF (LOWER(kwd) .EQ. mat_name(iPrp)) THEN
            mat%prp(iPrp) = u
            RETURN
         END IF
      END DO
      io%e = "Unable to find material property <"//TRIM(kwd)//">"

      RETURN
      END SUBROUTINE setMat
!#####################################################################
      FUNCTION FIND_MAT(mName) RESULT(s)
      CHARACTER(LEN=*), INTENT(IN) :: mName
      TYPE(matType), POINTER :: s

      INTEGER iMat

      DO iMat=1, SIZE(matLib)
         IF (matLib(iMat)%name .EQ. LOWER(mName)) THEN
            s => matLib(iMat)
            RETURN
         END IF
      END DO
      io%e = "Unable to find material <"//TRIM(mName)//">"

      RETURN
      END FUNCTION FIND_MAT
!---------------------------------------------------------------------      
      SUBROUTINE ADD_MAT(mat)
      TYPE(matType), INTENT(IN) :: mat

      INTEGER n
      TYPE(matType), ALLOCATABLE :: tMat(:)

      IF (ALLOCATED(matLib)) THEN
         n = SIZE(matLib)
         CALL MOVE_ALLOC(matLib, tMat)
         ALLOCATE(matLib(n+1))
         matLib(1:n) = tMat
         matLib(n+1) = mat
         DEALLOCATE(tMat)
      ELSE
         ALLOCATE(matLib(1))
         matLib = mat
      END IF

      RETURN
      END SUBROUTINE ADD_MAT
!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_MATMOD()

      TYPE(matType) mat
      TYPE(matType), POINTER :: mPtr

      mat = matType('Water',(/"Density","Conductivity","Specific heat",
     2   "Heat source","Dynamic viscosity"/),(/2D0,5D0,1D0,2D0,1D-1/)) 
      IF (mat%rho() .NE. 2D0) io%e = "Issue with newMat"
      io%o = "newMat: "//CLR("(PASSED)",3)
      CALL mat%set("Density",1D0)
      IF (mat%rho() .NE. 1D0) io%e = "Issue with mat%set"
      io%o = "mat%set: "//CLR("(PASSED)",3)
      IF (mat%E() .NE. 0D0) io%e = "Issue with mat%prp"
      io%o = "mat%prp: "//CLR("(PASSED)",3)
      CALL ADD_MAT(mat)
      io%o = "ADD_MAT: "//CLR("(PASSED)",3)
      mat = matType('Steel',(/"Density","Elasticity modulus","Diameter",
     2   "Poisson's ratio","Body force_z"/),(/1D2,1D3,1D-4,0.4D0,1D1/))
      CALL ADD_MAT(mat)
      mPtr => FIND_MAT('Water')
      IF (mPtr%rho() .NE. 1D0) io%e = "Issue with FIND_MAT"
      io%o = "FIND_MAT: "//CLR("(PASSED)",3)
     
      RETURN
      END SUBROUTINE TEST_MATMOD
      END MODULE MATMOD
