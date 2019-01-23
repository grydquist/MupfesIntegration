!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     For controlling the solution procesure. Whether to keep iterating
!     or go to the next time step. For usage and testing the class
!     see TEST_SCTMOD at the end of this file. 
!      
!--------------------------------------------------------------------
      MODULE SCTMOD
      USE LSTMOD
      IMPLICIT NONE

      TYPE sctType
!        Number of performed iterations
         INTEGER :: itr = 0
!        Maximum iteration for this eq.
         INTEGER :: maxItr = 10
!        Minimum iteration for this eq.
         INTEGER :: minItr = 1
!        dB reduction in residual
         REAL(KIND=8) :: dBr = -6D1
!        Residual at t=0
         REAL(KIND=8) :: zNorm = 0D0
!        First iteration residual
         REAL(KIND=8) pNorm
!        Current residual
         REAL(KIND=8) cNorm
!        Accepted relative tolerance
         REAL(KIND=8) :: tol = 1D64
      CONTAINS 
!        Whether the equation is satisfied
         PROCEDURE :: satisfied => satisfiedSct
!        Updates the norm enteries
         PROCEDURE :: upNorm => upNormSct
      END TYPE sctType

      INTERFACE sctType
         PROCEDURE :: newSct, newSctFl
      END INTERFACE 

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
      FUNCTION newSct(tol, dBr, minItr, maxItr) RESULT(sct)
      TYPE(sctType) :: sct
      REAL(KIND=8), INTENT(IN), OPTIONAL :: tol, dBr
      INTEGER, INTENT(IN), OPTIONAL :: minItr, maxItr

      IF (PRESENT(tol   )) sct%tol    = tol
      IF (PRESENT(dBr   )) sct%dBr    = dBr
      IF (PRESENT(minItr)) sct%minItr = minItr
      IF (PRESENT(maxItr)) sct%maxItr = maxItr
      
      RETURN
      END FUNCTION newSct
!---------------------------------------------------------------------
      FUNCTION newSctFl(lst) RESULT(sct)
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(sctType) :: sct

      TYPE(lstType), POINTER :: lPtr

      lPtr => lst%get(sct%minItr,"Min iterations",ll=1)
      lPtr => lst%get(sct%maxItr,"Max iterations",ll=1)
      lPtr => lst%get(sct%dBr,"Residual dB reduction",ul=0D0)
      lPtr => lst%get(sct%tol,"Tolerance",ll=0D0)

      RETURN
      END FUNCTION newSctFl
!---------------------------------------------------------------------
!     Updating the norm
      FUNCTION satisfiedSct(sct) RESULT(flag)
      CLASS(sctType), INTENT(IN) :: sct
      LOGICAL :: flag

      REAL(KIND=8) dBr

!     Exceeding maxItr overrides other conditions
      IF (sct%itr .GE. sct%maxItr) THEN
         flag = .TRUE. 
         RETURN
      END IF

!     Not satisfying minItr overrides other conditions
      IF (sct%itr .LT. sct%minItr) THEN
         flag = .FALSE. 
         RETURN
      END IF

      dBr = 2D1*LOG10(sct%cNorm/sct%pNorm)
      IF (dBr.LE.sct%dBr .AND. sct%cNorm.LE.sct%tol) THEN
         flag = .TRUE.
      ELSE
         flag = .FALSE.
      END IF

      RETURN
      END FUNCTION satisfiedSct
!---------------------------------------------------------------------
!     Updating the norm
      SUBROUTINE upNormSct(sct,nrm)
      CLASS(sctType), INTENT(INOUT) :: sct
      REAL(KIND=8), INTENT(IN) :: nrm

!     Increasing iteration by one
      sct%itr = sct%itr + 1

!     This should be exectued only once
      IF (ISZERO(sct%zNorm)) sct%zNorm = nrm

!     Keeping a non-dimensional version of nrm
      sct%cNorm = nrm/MAX(sct%zNorm,eps)

!     Setting pNorm to be the nrm associated with the first iteration
!     (in a non-dimensional form)
      IF (sct%itr .EQ. 1) sct%pNorm = sct%cNorm

      RETURN  
      END SUBROUTINE upNormSct
!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_SCTMOD(sct)
      TYPE(sctType), INTENT(OUT) :: sct

      sct  = sctType(tol=2D-2,dBr=-20D0)
      io%o = "newSct: "//CLR("(PASSED)",3)
      CALL sct%upNorm(1D3)
      IF (sct%zNorm.NE.1D3 .OR. sct%pNorm.NE.1D0) 
     2   io%e = "Issue with %upNorm"
      io%o = "sct%upNorm: "//CLR("(PASSED)",3)
      CALL sct%upNorm(1D2)
      IF (sct%satisfied()) io%e = "Issue with %satisfied"
      CALL sct%upNorm(1D1)
      IF (.NOT.sct%satisfied()) io%e = "Issue with %satisfied"
      io%o = "sct%satisifed: "//CLR("(PASSED)",3)
!     To reset 
      sct%itr   = 0 
      sct%zNorm = 0D0

      RETURN
      END SUBROUTINE TEST_SCTMOD
      END MODULE SCTMOD
