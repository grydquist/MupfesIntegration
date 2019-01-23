!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Contains element class. Represetation of element, shape function
!     calculation, Gauss's quadratures, shape function gradients. Also
!     for vtk types for IO. For usage and testing this module see
!     TEST_ELEMOD at the end of this file. 
!      
!--------------------------------------------------------------------
      MODULE ELEMOD
      USE IOMOD
      IMPLICIT NONE

!     Gauss points and their corresponding weights, upto 5 points
      REAL(KIND=8), PARAMETER :: gW(5,5)=RESHAPE((/2D0,0D0,0D0,0D0,0D0,
     2   1D0,1D0,0D0,0D0,0D0, 0.5555555555555556D0,0.8888888888888889D0,
     3   0.5555555555555556D0,0D0,0D0, 0.3478548451374538D0,
     4   0.6521451548625462D0,0.6521451548625462D0,0.3478548451374538D0,
     5   0D0, 0.236926885056189D0,0.4786286704993665D0,
     6   0.5688888888888889D0,0.4786286704993665D0,0.236926885056189D0/)
     7   ,(/5,5/))
      REAL(KIND=8), PARAMETER :: gXi(5,5)=RESHAPE((/0D0,0D0,0D0,0D0,0D0,
     2   -0.57735026918962584D0,0.57735026918962584D0,0D0,0D0,0D0, 
     3   -0.7745966692414834D0,0D0,0.7745966692414834D0,0D0,0D0,
     4   -0.86113631159405257D0,-0.33998104358485631D0,
     5   0.33998104358485631D0,0.86113631159405257D0,0D0,
     6   -0.90617984593866396D0,-0.53846931010568311D0,0D0,
     7   0.53846931010568311D0,0.90617984593866396D0/),(/5,5/))

!     Types of accepted elements
!     Linear (1D), triangle (2D), tetrahedral (3D), bilinear (2D), quad
!     (1D), biquad (2D), brick (3D), general NURBS (1-3D)
      INTEGER, PARAMETER :: eType_NA = 100, eType_LIN = 101, 
     2   eType_TRI = 102, eType_TET = 103, eType_BIL = 104, 
     3   eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107, 
     4   eType_NRB = 108, eType_WDG = 109

      TYPE eleType
!        Whether the shape function is linear
         LOGICAL :: lShpF
!        Element type
         INTEGER :: eType = eType_NA
!        Number of nodes (control points) in a single element
         INTEGER :: eNoN 
!        Number of element face
         INTEGER :: nEf
!        Number of Gauss points for integration
         INTEGER :: nG
!        the element type recognized by VTK format
         INTEGER :: vtkType
!        Gauss weights: nG
         REAL(KIND=8), ALLOCATABLE :: w(:)
!        Local location of Gauss points: nsd x nG
         REAL(KIND=8), ALLOCATABLE :: xi(:,:)
!        Parent shape function at Gauss points: eNoN x nG
         REAL(KIND=8), ALLOCATABLE :: N(:,:)
!        Parent shape functions gradient at Gauss pnts: nsd x eNoN x nG
         REAL(KIND=8), ALLOCATABLE :: Nx(:,:,:)
      CONTAINS
!        Computes shape functions gradient at a Gauss point
         PROCEDURE :: dNdx => dNdxEle
!        Computes shape function at a given point
         PROCEDURE :: NAtx => NAtxEle
!        To deallocate an element
         FINAL :: freeEle

!        Assigns all shape functions, w, xi, N, and Nx
         PROCEDURE, PRIVATE :: set => setEle
!        Computes shape function at a given local point
         PROCEDURE, PRIVATE :: NAtxi => NAtxiEle
!        Computes shape function gradient at a given local point
         PROCEDURE, PRIVATE :: NxAtxi => NxAtxiEle
      END TYPE eleType

      INTERFACE eleType
         PROCEDURE :: newEle
      END INTERFACE eleType
!---------------------------------------------------------------------
!     Number of spatial dimensions
      INTEGER :: nsd = 3

      CONTAINS
!####################################################################

!     IMPLEMENTATION

!####################################################################
!     Finds the element type given element number of nodes as well as
!     the dimension of space that this element/facet belongs to
      FUNCTION ELEMENT_TYPE(iNSD, eNoN) RESULT(eType)
      INTEGER, INTENT(IN) :: iNSD, eNoN
      INTEGER eType

      eType = eType_NA
      IF (iNSD .EQ. 3) THEN
         SELECT CASE (eNoN)
         CASE(8)
            eType = eType_BRK
         CASE(6)
            eType = eType_WDG
         CASE(4)
            eType = eType_TET
         CASE DEFAULT
           io%e = "Unable to identify combination of nsd and eNoN"
         END SELECT
      ELSE IF (iNSD .EQ. 2) THEN
         SELECT CASE (eNoN)
         CASE(3)
            eType = eType_TRI
         CASE(4)
            eType = eType_BIL
         CASE(9)
            eType = eType_BIQ
         CASE DEFAULT
            io%e = "Unable to identify combination of nsd and eNoN"
         END SELECT
      ELSE IF (iNSD .EQ. 1) THEN
         SELECT CASE (eNoN)
         CASE(2)
            eType = eType_LIN
         CASE(3)
            eType = eType_QUD
         CASE DEFAULT
            io%e = "Unable to identify combination of nsd and eNoN"
         END SELECT
      END IF

      RETURN
      END FUNCTION ELEMENT_TYPE

!####################################################################

      FUNCTION newEle(eType, eNoN, nG) RESULT(e)
      TYPE(eleType) :: e
      INTEGER, INTENT(IN) :: eType
      INTEGER, INTENT(IN), OPTIONAL :: eNoN, nG

      e%eType = eType
      SELECT CASE (eType)
      CASE (eType_BRK)
         e%eNoN    = 8
         e%nG      = 8
         e%vtkType = 12
         e%nEf     = 6
         e%lShpF   = .FALSE.
      CASE (eType_WDG)
         e%eNoN    = 6
         e%nG      = 6
         e%vtkType = 13
         e%nEf     = 3
         e%lShpF   = .FALSE.
      CASE (eType_TET)
         e%eNoN    = 4
         e%nG      = 4
         e%vtkType = 10
         e%nEf     = 4
         e%lShpF   = .TRUE.
      CASE (eType_TRI)
         e%eNoN    = 3
         e%nG      = 3
         e%vtkType = 5
         e%nEf     = 3
         e%lShpF   = .TRUE.
      CASE (eType_BIL)
         e%eNoN    = 4
         e%nG      = 4
         e%vtkType = 9
         e%nEf     = 4
         e%lShpF   = .FALSE.
      CASE (eType_BIQ)
         e%eNoN    = 9
         e%nG      = 9
         e%vtkType = 28
         e%nEf     = 4
         e%lShpF   = .FALSE.
      CASE(eType_LIN)
         e%eNoN    = 2
         e%nG      = 2
         e%vtkType = 0
         e%nEf     = 2
         e%lShpF   = .TRUE.
      CASE(eType_QUD)
         e%eNoN    = 3
         e%nG      = 3
         e%vtkType = 0
         e%nEf     = 2
         e%lShpF   = .FALSE.
      CASE (eType_NRB)
         IF (.NOT.PRESENT(nG) .OR. .NOT.PRESENT(eNoN)) 
     2      io%e = "newEle: eNoN and nG must be specified for NURBS"

         e%eNoN  = eNoN
         e%nG    = nG
         e%lShpF = .FALSE.
         IF (nsd .EQ. 3) THEN
            e%vtkType = 12
            e%nEf     = 6
         ELSE IF (nsd .EQ. 2) THEN
            e%vtkType = 9
            e%nEf     = 4
         ELSE
            io%e = "newEle: Undefined nsd"
         END IF
      CASE DEFAULT
         io%e = "Unable to identify eType in newEle"
      END SELECT
      
      CALL e%set()

      RETURN
      END FUNCTION newEle
!---------------------------------------------------------------------
!     Assigning values to w, xi, N and Nx
      PURE SUBROUTINE setEle(e)
      CLASS(eleType), INTENT(INOUT) :: e

      INTEGER g
      REAL(KIND=8) s, t, uz, lz

      ALLOCATE(e%w(e%nG), e%xi(nsd,e%nG), e%N(e%eNoN,e%nG), 
     2   e%Nx(nsd,e%eNoN,e%nG))

      IF (e%eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(e%eType)
      CASE(eType_BRK)
         e%w = 1D0
         s   =  1D0/SQRT(3D0)
         t   = -1D0/SQRT(3D0)
         e%xi(1,1) = s; e%xi(2,1) = s; e%xi(3,1) = s
         e%xi(1,2) = t; e%xi(2,2) = s; e%xi(3,2) = s
         e%xi(1,3) = t; e%xi(2,3) = s; e%xi(3,3) = t
         e%xi(1,4) = s; e%xi(2,4) = s; e%xi(3,4) = t
         e%xi(1,5) = s; e%xi(2,5) = t; e%xi(3,5) = s
         e%xi(1,6) = t; e%xi(2,6) = t; e%xi(3,6) = s
         e%xi(1,7) = t; e%xi(2,7) = t; e%xi(3,7) = t
         e%xi(1,8) = s; e%xi(2,8) = t; e%xi(3,8) = t
      CASE(eType_TET)
         e%w = 1D0/24D0
         s   = (5D0 + 3D0*SQRT(5D0))/2D1
         t   = (5D0 -     SQRT(5D0))/2D1
         e%xi(1,1) = s; e%xi(2,1) = t; e%xi(3,1) = t
         e%xi(1,2) = t; e%xi(2,2) = s; e%xi(3,2) = t
         e%xi(1,3) = t; e%xi(2,3) = t; e%xi(3,3) = s
         e%xi(1,4) = t; e%xi(2,4) = t; e%xi(3,4) = t
      CASE(eType_WDG)
         e%w =  1D0/6D0
         s   =  2D0/3D0
         t   =  1D0/6D0
         uz  =  1D0/SQRT(3D0)
         lz  = -1D0/SQRT(3D0)
         e%xi(1,1) = s; e%xi(2,1) = t; e%xi(3,1) = lz
         e%xi(1,2) = t; e%xi(2,2) = s; e%xi(3,2) = lz
         e%xi(1,3) = t; e%xi(2,3) = t; e%xi(3,3) = lz
         e%xi(1,4) = s; e%xi(2,4) = t; e%xi(3,4) = uz
         e%xi(1,5) = t; e%xi(2,5) = s; e%xi(3,5) = uz
         e%xi(1,6) = t; e%xi(2,6) = t; e%xi(3,6) = uz
!     2D elements         
      CASE(eType_TRI)
         e%w = 1D0/6D0
         s   = 2D0/3D0
         t   = 1D0/6D0
         e%xi(1,1) = s; e%xi(2,1) = t
         e%xi(1,2) = t; e%xi(2,2) = s
         e%xi(1,3) = t; e%xi(2,3) = t
      CASE(eType_BIL)
         e%w = 1D0
         s   =  1D0/SQRT(3D0)
         t   = -1D0/SQRT(3D0)
         e%xi(1,1) = s; e%xi(2,1) = s
         e%xi(1,2) = t; e%xi(2,2) = s
         e%xi(1,3) = t; e%xi(2,3) = t
         e%xi(1,4) = s; e%xi(2,4) = t
      CASE(eType_BIQ)
         e%w(1) = 25D0/81D0
         e%w(2) = 25D0/81D0
         e%w(3) = 25D0/81D0
         e%w(4) = 25D0/81D0
         e%w(5) = 40D0/81D0
         e%w(6) = 40D0/81D0
         e%w(7) = 40D0/81D0
         e%w(8) = 40D0/81D0
         e%w(9) = 64D0/81D0
         s      = SQRT(6D-1)
         e%xi(1,1) =  -s; e%xi(2,1) =  -s
         e%xi(1,2) =   s; e%xi(2,2) =  -s
         e%xi(1,3) =   s; e%xi(2,3) =   s
         e%xi(1,4) =  -s; e%xi(2,4) =   s
         e%xi(1,5) = 0D0; e%xi(2,5) =  -s
         e%xi(1,6) =   s; e%xi(2,6) = 0D0
         e%xi(1,7) = 0D0; e%xi(2,7) =   s
         e%xi(1,8) =  -s; e%xi(2,8) = 0D0
         e%xi(1,9) = 0D0; e%xi(2,9) = 0D0
!     1D elements         
      CASE(eType_LIN)
         e%w = 1D0
         s   = 1D0/SQRT(3D0)
         e%xi(1,1) = -s
         e%xi(1,2) =  s
      CASE(eType_QUD)
         e%w(1) = 5D0/9D0
         e%w(2) = 5D0/9D0
         e%w(3) = 8D0/9D0
         s = SQRT(6D-1)
         e%xi(1,1) = -s
         e%xi(1,2) =  s
         e%xi(1,3) = 0D0
      END SELECT
      
      DO g=1, e%nG
         e%N(:,g)    = e%NAtxi (e%xi(:,g))
         e%Nx(:,:,g) = e%NxAtxi(e%xi(:,g))
      END DO

      RETURN
      END SUBROUTINE setEle
!---------------------------------------------------------------------
!     Computing d(N_a)/d(x) at a Gauss point
      PURE SUBROUTINE dNdxEle(e, g, x, Nx, Jac, ks)
      CLASS(eleType), INTENT(IN) :: e
      INTEGER, INTENT(IN) :: g
      REAL(KIND=8), INTENT(IN) :: x(nsd,e%eNoN)
      REAL(KIND=8), INTENT(OUT) :: Nx(nsd,e%eNoN), Jac
      REAL(KIND=8), INTENT(OUT), OPTIONAL :: ks(nsd,nsd)

      INTEGER a
      REAL(KIND=8) xXi(nsd,nsd), xiX(nsd,nsd), iJac

!     The idea is as follows: if we express x as it is a variable
!     defined in the parent coordinate system and discretized using 
!     our shape functions, then x(xi) = \sum_a N_a(xi) x_a. Taking
!     derivative of this relationship gives: d(x)/d(xi) = xXi = \sum_a
!     d(N_a)/d(xi) x_a. Since d(N_a)/d(xi) is available from element
!     class, xXi can be readily computed. Also, we have: d(xi)/d(x) =
!     xiX = (d(x)/d(xi))^-1 which is computed at the next step. Finally
!     to computed Nx = d(N_a)/d(x) we use chain rule: d(N_a)/d(x) =
!     d(N_a)/d(xi) x d(xi)/d(x). 

      Nx  = 0D0
      xXi = 0D0
      IF (nsd .EQ. 2) THEN
         DO a=1, e%eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*e%Nx(1,a,g)
            xXi(:,2) = xXi(:,2) + x(:,a)*e%Nx(2,a,g)
         END DO

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)
         iJac = 1D0/Jac

         xiX(1,1) =  xXi(2,2)*iJac
         xiX(1,2) = -xXi(1,2)*iJac
         xiX(2,1) = -xXi(2,1)*iJac
         xiX(2,2) =  xXi(1,1)*iJac

         DO a=1, e%eNoN
            Nx(1,a) = Nx(1,a) + e%Nx(1,a,g)*xiX(1,1) 
     2                        + e%Nx(2,a,g)*xiX(2,1)
            Nx(2,a) = Nx(2,a) + e%Nx(1,a,g)*xiX(1,2) 
     2                        + e%Nx(2,a,g)*xiX(2,2)
         END DO

         IF (PRESENT(ks)) THEN
            ks(1,1) = xiX(1,1)*xiX(1,1) + xiX(2,1)*xiX(2,1)
            ks(1,2) = xiX(1,1)*xiX(1,2) + xiX(2,1)*xiX(2,2)
            ks(2,2) = xiX(1,2)*xiX(1,2) + xiX(2,2)*xiX(2,2)
            ks(2,1) = ks(1,2)
         END IF
      ELSE
         DO a=1, e%eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*e%Nx(1,a,g)
            xXi(:,2) = xXi(:,2) + x(:,a)*e%Nx(2,a,g)
            xXi(:,3) = xXi(:,3) + x(:,a)*e%Nx(3,a,g)
         END DO
         
         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2       + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3       + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4       - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5       - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6       - xXi(1,3)*xXi(2,2)*xXi(3,1)
         iJac = 1D0/Jac

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))*iJac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))*iJac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))*iJac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))*iJac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))*iJac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))*iJac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))*iJac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))*iJac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))*iJac
        
         DO a=1, e%eNoN
            Nx(1,a) = Nx(1,a) + e%Nx(1,a,g)*xiX(1,1) 
     2                        + e%Nx(2,a,g)*xiX(2,1) 
     3                        + e%Nx(3,a,g)*xiX(3,1)
            
            Nx(2,a) = Nx(2,a) + e%Nx(1,a,g)*xiX(1,2) 
     2                        + e%Nx(2,a,g)*xiX(2,2) 
     3                        + e%Nx(3,a,g)*xiX(3,2)
            
            Nx(3,a) = Nx(3,a) + e%Nx(1,a,g)*xiX(1,3) 
     2                        + e%Nx(2,a,g)*xiX(2,3) 
     3                        + e%Nx(3,a,g)*xiX(3,3)
         END DO
         
         IF (PRESENT(ks)) THEN
            ks(1,1) = xiX(1,1)*xiX(1,1) 
     2              + xiX(2,1)*xiX(2,1) 
     3              + xiX(3,1)*xiX(3,1)
            ks(1,2) = xiX(1,2)*xiX(1,1) 
     2              + xiX(2,2)*xiX(2,1)
     3              + xiX(3,2)*xiX(3,1)
            ks(1,3) = xiX(1,3)*xiX(1,1) 
     2              + xiX(2,3)*xiX(2,1)
     3              + xiX(3,3)*xiX(3,1)
            ks(2,2) = xiX(1,2)*xiX(1,2) 
     2              + xiX(2,2)*xiX(2,2)
     3              + xiX(3,2)*xiX(3,2)
            ks(2,3) = xiX(1,2)*xiX(1,3) 
     2              + xiX(2,2)*xiX(2,3)
     3              + xiX(3,2)*xiX(3,3)
            ks(3,3) = xiX(1,3)*xiX(1,3) 
     2              + xiX(2,3)*xiX(2,3)
     3              + xiX(3,3)*xiX(3,3)
            ks(2,1) = ks(1,2)
            ks(3,1) = ks(1,3)
            ks(3,2) = ks(2,3)
         END IF
      END IF

      RETURN
      END SUBROUTINE dNdxEle
!--------------------------------------------------------------------
!     Computes shape function N of an element for an arbitrary point:
!     x. If the point belongs to the element, true value is returned,
!     otherwise false is returned. N remains unchanged in case
!     of false return
      PURE FUNCTION NAtxEle(e, x, xl) RESULT(N)
      CLASS(eleType), INTENT(IN) :: e
      REAL(KIND=8), INTENT(IN) :: x(nsd), xl(nsd,e%eNoN)
      REAL(KIND=8) :: N(e%eNoN)

      INTEGER ia
      REAL(KIND=8) xXi(nsd,nsd), A(nsd,nsd), xr(nsd), Jac, iJac,
     2   xi(nsd)

!     The idea is as follows: assuming xi is a linear transformation of
!     x then: xi = A*x + B with A nsd x nsd matrix and B nsd x 1 vector.
!     Taking derivative of this relationship, A = d(xi)/d(x) which can
!     be computed using a procedure described above. To compute B, we
!     resort to a given known point, i.e. the first Gauss point at which 
!     we have x_1 = sum_a N_a(xi_1) x_a and xi_1. Thus B = xi_1 -
!     A*x_1. This relation in combination with xi = A*x + B yields: 
!     xi = A*(x - x_1) + xi_1 = A*xr + xi_1, which is implemented in 
!     practice. Given xi, then shape functions are computed by a call 
!     to NAtXi. 

      xXi = 0D0
      xr  = x 
      DO ia=1, e%eNoN
         xXi(:,1) = xXi(:,1) + xl(:,ia)*e%Nx(1,ia,1)
         xXi(:,2) = xXi(:,2) + xl(:,ia)*e%Nx(2,ia,1)
         xXi(:,3) = xXi(:,3) + xl(:,ia)*e%Nx(3,ia,1)
         xr       = xr       - xl(:,ia)*e%N (  ia,1)
      END DO

      Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2    + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3    + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4    - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5    - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6    - xXi(1,3)*xXi(2,2)*xXi(3,1)
      iJac = 1D0/Jac

      A(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))*iJac
      A(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))*iJac
      A(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))*iJac
      A(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))*iJac
      A(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))*iJac
      A(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))*iJac
      A(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))*iJac
      A(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))*iJac
      A(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))*iJac


      xi(1) = A(1,1)*xr(1) + A(1,2)*xr(2) + A(1,3)*xr(3) + e%xi(1,1)
      xi(2) = A(2,1)*xr(1) + A(2,2)*xr(2) + A(2,3)*xr(3) + e%xi(2,1)
      xi(3) = A(3,1)*xr(1) + A(3,2)*xr(2) + A(3,3)*xr(3) + e%xi(3,1)
 
      N = e%NAtxi(xi)

      RETURN
      END FUNCTION NAtxEle
!--------------------------------------------------------------------
!     Computes shape functions N at a given local point xi.
      PURE FUNCTION NAtxiEle(e, xi) RESULT(N)
      CLASS(eleType), INTENT(IN) :: e
      REAL(KIND=8), INTENT(IN) :: xi(nsd)
      REAL(KIND=8) :: N(e%eNoN)

      REAL(KIND=8) s, t, mx, my, ux, uy, uz, lx, ly, lz
      
      IF (e%eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(e%eType)
      CASE(eType_BRK)
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         uz = 1D0 + xi(3)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         lz = 1D0 - xi(3)

         N(1) = ux*uy*uz/8D0
         N(2) = lx*uy*uz/8D0
         N(3) = lx*uy*lz/8D0
         N(4) = ux*uy*lz/8D0
         N(5) = ux*ly*uz/8D0
         N(6) = lx*ly*uz/8D0
         N(7) = lx*ly*lz/8D0
         N(8) = ux*ly*lz/8D0
      CASE(eType_TET)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1D0 - xi(1) - xi(2) - xi(3)
      CASE(eType_WDG)
         ux = xi(1)
         uy = xi(2)
         uz = 1D0 - xi(1) - xi(2)
         s  = (1D0 + xi(3))/2D0
         t  = (1D0 - xi(3))/2D0
         
         N(1) = ux*t
         N(2) = uy*t
         N(3) = uz*t
         N(4) = ux*s
         N(5) = uy*s
         N(6) = uz*s

!     2D elements         
      CASE(eType_TRI)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = 1D0 - xi(1) - xi(2)
      CASE(eType_BIL)
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         
         N(1) = ux*uy/4D0
         N(2) = lx*uy/4D0
         N(3) = lx*ly/4D0
         N(4) = ux*ly/4D0
      CASE(eType_BIQ)
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         mx = xi(1)
         my = xi(2)
         
         N(1) =  mx*lx*my*ly/4D0
         N(2) = -mx*ux*my*ly/4D0
         N(3) =  mx*ux*my*uy/4D0
         N(4) = -mx*lx*my*uy/4D0
         N(5) = -lx*ux*my*ly/2D0
         N(6) =  mx*ux*ly*uy/2D0
         N(7) =  lx*ux*my*uy/2D0
         N(8) = -mx*lx*ly*uy/2D0
         N(9) =  lx*ux*ly*uy

!     1D elements         
      CASE(eType_LIN)
         N(1) = (1D0 - xi(1))/2D0
         N(2) = (1D0 + xi(1))/2D0
      CASE(eType_QUD)
         N(1) = -xi(1)*(1D0 - xi(1))/2D0
         N(2) =  xi(1)*(1D0 + xi(1))/2D0
         N(3) = (1D0 - xi(1))*(1D0 + xi(1))
      END SELECT

      RETURN
      END FUNCTION NAtxiEle
!--------------------------------------------------------------------
!     Computes shape functions gradient Nx at a given local point xi.
      PURE FUNCTION NxAtxiEle(e, xi) RESULT(Nx)
      CLASS(eleType), INTENT(IN) :: e
      REAL(KIND=8), INTENT(IN) :: xi(nsd)
      REAL(KIND=8) :: Nx(nsd,e%eNoN)

      REAL(KIND=8) s, t, mx, my, ux, uy, uz, lx, ly, lz
      
      IF (e%eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(e%eType)
      CASE(eType_BRK)
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         uz = 1D0 + xi(3)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         lz = 1D0 - xi(3)

         Nx(1,1) =  uy*uz/8D0
         Nx(2,1) =  ux*uz/8D0
         Nx(3,1) =  ux*uy/8D0
         Nx(1,2) = -uy*uz/8D0
         Nx(2,2) =  lx*uz/8D0
         Nx(3,2) =  lx*uy/8D0
         Nx(1,3) = -uy*lz/8D0
         Nx(2,3) =  lx*lz/8D0
         Nx(3,3) = -lx*uy/8D0
         Nx(1,4) =  uy*lz/8D0
         Nx(2,4) =  ux*lz/8D0
         Nx(3,4) = -ux*uy/8D0
         Nx(1,5) =  ly*uz/8D0
         Nx(2,5) = -ux*uz/8D0
         Nx(3,5) =  ux*ly/8D0
         Nx(1,6) = -ly*uz/8D0
         Nx(2,6) = -lx*uz/8D0
         Nx(3,6) =  lx*ly/8D0
         Nx(1,7) = -ly*lz/8D0
         Nx(2,7) = -lx*lz/8D0
         Nx(3,7) = -lx*ly/8D0
         Nx(1,8) =  ly*lz/8D0
         Nx(2,8) = -ux*lz/8D0
         Nx(3,8) = -ux*ly/8D0
      CASE(eType_TET)
         Nx(1,1) =  1D0
         Nx(2,1) =  0D0
         Nx(3,1) =  0D0
         Nx(1,2) =  0D0
         Nx(2,2) =  1D0
         Nx(3,2) =  0D0
         Nx(1,3) =  0D0
         Nx(2,3) =  0D0
         Nx(3,3) =  1D0
         Nx(1,4) = -1D0
         Nx(2,4) = -1D0
         Nx(3,4) = -1D0
      CASE(eType_WDG)
         ux = xi(1)
         uy = xi(2)
         uz = 1D0 - xi(1) - xi(2)
         s  = (1D0 + xi(3))/2D0
         t  = (1D0 - xi(3))/2D0
         
         Nx(1,1) =  t
         Nx(2,1) =  0D0
         Nx(3,1) = -ux/2D0
         Nx(1,2) =  0D0
         Nx(2,2) =  t
         Nx(3,2) = -uy/2D0
         Nx(1,3) = -t
         Nx(2,3) = -t
         Nx(3,3) = -uz/2D0
         Nx(1,4) =  s
         Nx(2,4) =  0D0
         Nx(3,4) =  ux/2D0
         Nx(1,5) =  0D0
         Nx(2,5) =  s
         Nx(3,5) =  uy/2D0
         Nx(1,6) = -s
         Nx(2,6) = -s
         Nx(3,6) =  uz/2D0

!     2D elements         
      CASE(eType_TRI)
         Nx(1,1) =  1D0
         Nx(2,1) =  0D0
         Nx(1,2) =  0D0
         Nx(2,2) =  1D0
         Nx(1,3) = -1D0
         Nx(2,3) = -1D0
      CASE(eType_BIL)
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         
         Nx(1,1) =  uy/4D0
         Nx(2,1) =  ux/4D0
         Nx(1,2) = -uy/4D0
         Nx(2,2) =  lx/4D0
         Nx(1,3) = -ly/4D0
         Nx(2,3) = -lx/4D0
         Nx(1,4) =  ly/4D0
         Nx(2,4) = -ux/4D0
      CASE(eType_BIQ)
         ux = 1D0 + xi(1)
         uy = 1D0 + xi(2)
         lx = 1D0 - xi(1)
         ly = 1D0 - xi(2)
         mx = xi(1)
         my = xi(2)
         
         Nx(1,1) =  (lx - mx)*my*ly/4D0
         Nx(2,1) =  (ly - my)*mx*lx/4D0
         Nx(1,2) = -(ux + mx)*my*ly/4D0
         Nx(2,2) = -(ly - my)*mx*ux/4D0
         Nx(1,3) =  (ux + mx)*my*uy/4D0
         Nx(2,3) =  (uy + my)*mx*ux/4D0
         Nx(1,4) = -(lx - mx)*my*uy/4D0
         Nx(2,4) = -(uy + my)*mx*lx/4D0
         Nx(1,5) = -(lx - ux)*my*ly/2D0
         Nx(2,5) = -(ly - my)*lx*ux/2D0
         Nx(1,6) =  (ux + mx)*ly*uy/2D0
         Nx(2,6) =  (ly - uy)*mx*ux/2D0
         Nx(1,7) =  (lx - ux)*my*uy/2D0
         Nx(2,7) =  (uy + my)*lx*ux/2D0
         Nx(1,8) = -(lx - mx)*ly*uy/2D0
         Nx(2,8) = -(ly - uy)*mx*lx/2D0
         Nx(1,9) =  (lx - ux)*ly*uy
         Nx(2,9) =  (ly - uy)*lx*ux

!     1D elements         
      CASE(eType_LIN)
         Nx(1,1) = -5D-1
         Nx(1,2) =  5D-1
      CASE(eType_QUD)
         Nx(1,1) = -5D-1 + xi(1)
         Nx(1,2) =  5D-1 + xi(1)
         Nx(1,3) = -2D0*xi(1)
      END SELECT

      RETURN
      END FUNCTION NxAtxiEle
!--------------------------------------------------------------------
      SUBROUTINE freeEle(this)
      TYPE(eleType) :: this

      this%eType = eType_NA
      IF (ALLOCATED(this%Nx)) DEALLOCATE(this%Nx)
      IF (ALLOCATED(this%xi)) DEALLOCATE(this%xi)
      IF (ALLOCATED(this%w))  DEALLOCATE(this%w)
      IF (ALLOCATED(this%N))  DEALLOCATE(this%N)
      
      RETURN 
      END SUBROUTINE freeEle

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_ELEMOD()

      INTEGER :: g
      REAL(KIND=8) :: x(3,4), xp(3), Nx(3,4), Jac, N(4)
      TYPE(eleType) e

      e = eleType(eType_TET)
      IF (.NOT.e%lShpF .OR. e%eNoN.NE.4 .OR. e%nEf.NE.4 .OR. e%nG.NE.4
     2   .OR. e%vtkType.NE.10) io%e = "Issue with newEle[static]"
      io%o = "newEle (static): "//CLR("(PASSED)",3)
      IF (ANY(e%w.NE.1D0/24D0) .OR. SUM(e%xi(1,:)).NE.1D0 .OR.
     2   SUM(e%N(1,:)).NE.1D0 .OR. ANY(SUM(e%Nx(:,:,1),2).NE.0D0)) 
     3   io%e = "Issue with newEle[dynamic]"
      io%o = "newEle (dynamic): "//CLR("(PASSED)",3)
      g = 1
      x = RESHAPE((/1D0,0D0,0D0,0D0,2D0,0D0,0D0,0D0,3D0,0D0,0D0,0D0/),
     2   (/3,4/))
      CALL e%dNdx(g, x, Nx, Jac)
      IF (ANY(Nx(:,2).NE.(/0D0,0.5D0,0D0/)) .OR. Jac.NE.6D0) io%e = 
     2   "Issue with ele%dNdx(eType_TET)"
      io%o = "ele%dNdx: "//CLR("(PASSED)",3)
      xp = (/0.25,0.25,0.25/)
      N  = e%NAtx(xp,x)
      IF (.NOT.ISZERO(N,(/0.25D0,0.125D0,0.25D0/3D0,1D0-1.375D0/3D0/)))
     2    io%e = "Issue with ele%NAtx(eType_TET)"
      io%o = "ele%NAtx: "//CLR("(PASSED)",3)

      RETURN
      END SUBROUTINE TEST_ELEMOD
      END MODULE ELEMOD
