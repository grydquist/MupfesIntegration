!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     This is for solving advection-diffusion-reaction (adr) equation.
!     For usage and testing see TEST_ADRMOD at the end of this file.
!      
!---------------------------------------------------------------------
      MODULE ADRMOD
      USE EQMOD
      IMPLICIT NONE

!     Is used to initialize adr equation
      TYPE(eqSpType), PARAMETER, PRIVATE :: eqSp = eqSpType(
     1   nVar = 1,
     2   itgO = 1,
     3   lsAl = LS_TYPE_GMRES,
     4   sym  = "AD")

!     Equation type. Extending solution control 
      TYPE, EXTENDS(eqType) :: adrType
         PRIVATE
!        Temprature
         TYPE(gVarType), POINTER :: T => NULL()
!        Velocity, from NS
         TYPE(varType), POINTER :: U => NULL()
      CONTAINS
!        Setups all structure
         PROCEDURE :: setup => setupAdr
!        Evaulate the governing equation
         PROCEDURE :: eval3 => eval3Adr
         PROCEDURE :: eval2 => eval2Adr
!        Evaulate the governing equation boundary contribution
         PROCEDURE :: bEval => bEvalAdr
      END TYPE adrType

      INTERFACE adrType
         PROCEDURE :: newAdr, newAdrFl
      END INTERFACE adrType

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
!     This is for solving adr equation in a fluid or solid.
      FUNCTION newAdr(dmn, sct, mName, nBc, ls, U) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: mName
      INTEGER, INTENT(IN) :: nBc
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      CLASS(varType), INTENT(IN), OPTIONAL, TARGET :: U
      TYPE(adrType) :: eq
      
      CALL eq%new(eqSp, dmn, sct, mName, nBc, ls)
      IF (PRESENT(U)) eq%U => U

      RETURN
      END FUNCTION newAdr
!---------------------------------------------------------------------
      FUNCTION newAdrFl(dmn, lst, U) RESULT(eq)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst
      CLASS(varType), INTENT(IN), OPTIONAL, TARGET :: U
      TYPE(adrType) :: eq
      
      CALL eq%new(eqSp, dmn, lst)
      IF (PRESENT(U)) eq%U => U

      RETURN
      END FUNCTION newAdrFl
!---------------------------------------------------------------------
      SUBROUTINE setupAdr(eq, var)
      CLASS(adrType), INTENT(INOUT), TARGET :: eq
      TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
      
      IF (PRESENT(var)) THEN
         IF (SIZE(var) .NE. 1) io%e = "setupAdr: Invalid var size"
         IF (var(1)%dof .NE. 1) io%e = "setupAdr: Invalid var nU"
         eq%var => var
      ELSE
         eq%var(1) = gVarType(1,'Temperature',eq%dmn)
      END IF
      eq%T => eq%var(1)

      IF (eq%mat%rho() .LE. 0D0) io%e = "setupAdr: rho <= 0"
      IF (eq%mat%Cp() .LE. 0D0)  io%e = "setupAdr: Cp <= 0"
      IF (eq%mat%K() .LE. 0D0)   io%e = "setupAdr: K <= 0"

      RETURN
      END SUBROUTINE setupAdr
!---------------------------------------------------------------------
      PURE SUBROUTINE eval3Adr(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(adrType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
      
      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,3D0,1D0/)
      INTEGER a, b, Ac
      REAL(KIND=8) T1, amd, wl, wr, Td, Tx(nsd), tauM, kU, 
     2   kssk, nTx, Tp, udTx, udNx(eNoN), u(nsd), C, alp, q

      T1  = eq%itg%af*eq%itg%gam*dt
      amd = eq%itg%am/T1
      C   = eq%mat%Cp()*eq%mat%rho()
      alp = eq%mat%K()/C
      q   = eq%mat%q()/C
      wr  = w*C
      wl  = wr*T1

      Td  = 0D0
      Tx  = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         Td = Td + N(a)*eq%T%A%s(Ac)
      
         Tx(1) = Tx(1) + Nx(1,a)*eq%T%s(Ac)
         Tx(2) = Tx(2) + Nx(2,a)*eq%T%s(Ac)
         Tx(3) = Tx(3) + Nx(3,a)*eq%T%s(Ac)
      END DO

!     This is adr equation in a fluid
      IF (ASSOCIATED(eq%U)) THEN
         u = 0D0
         DO a=1, eNoN
            Ac   = eqN(a)
            u(1) = u(1) + N(a)*eq%U%v(1,Ac)
            u(2) = u(2) + N(a)*eq%U%v(2,Ac)
            u(3) = u(3) + N(a)*eq%U%v(3,Ac)
         END DO
     
         IF (ASSOCIATED(eq%dmn%Um)) THEN
            DO a=1, eNoN
               Ac = eqN(a)
               u(1) = u(1) - N(a)*eq%dmn%Um%v(1,Ac)
               u(2) = u(2) - N(a)*eq%dmn%Um%v(2,Ac)
               u(3) = u(3) - N(a)*eq%dmn%Um%v(3,Ac)
            END DO
         END IF

         kU = u(1)*u(1)*ks(1,1) + u(2)*u(1)*ks(2,1) + u(3)*u(1)*ks(3,1) 
     2      + u(1)*u(2)*ks(1,2) + u(2)*u(2)*ks(2,2) + u(3)*u(2)*ks(3,2) 
     3      + u(1)*u(3)*ks(1,3) + u(2)*u(3)*ks(2,3) + u(3)*u(3)*ks(3,3)

         kssk = ks(1,1)*ks(1,1) + ks(2,1)*ks(2,1) + ks(3,1)*ks(3,1) 
     2        + ks(1,2)*ks(1,2) + ks(2,2)*ks(2,2) + ks(3,2)*ks(3,2)
     3        + ks(1,3)*ks(1,3) + ks(2,3)*ks(2,3) + ks(3,3)*ks(3,3)
         
         nTx = ks(1,1)*Tx(1)*Tx(1) 
     2       + ks(2,2)*Tx(2)*Tx(2) 
     3       + ks(3,3)*Tx(3)*Tx(3) 
     4       + (ks(1,2) + ks(2,1))*Tx(1)*Tx(2)
     5       + (ks(1,3) + ks(3,1))*Tx(1)*Tx(3)
     6       + (ks(2,3) + ks(3,2))*Tx(2)*Tx(3)
         IF (ISZERO(nTx)) nTx = eps
      
         udTx = u(1)*Tx(1)   + u(2)*Tx(2)   + u(3)*Tx(3)
         udNx = u(1)*Nx(1,:) + u(2)*Nx(2,:) + u(3)*Nx(3,:)
         Tp   = ABS(Td + udTx)
!     Total viscosity is the physical + discontnuity capturing viscosity
         alp  = alp + Tp/SQRT(nTx)/2D0
         tauM = ct(4)/SQRT(ct(1)/dt/dt + ct(2)*kU + ct(3)*alp*alp*kssk)
         Tp   = -tauM*(Td + udTx)
      ELSE
         udTx = 0D0
         udNx = 0D0
         tauM = 0D0
         Tp   = 0D0
      END IF

!     Here the loop is started for constructing left and right hand side
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + wr*(N(a)*(Td + udTx) - udNx(a)*Tp
     2      + alp*(Nx(1,a)*Tx(1)+Nx(2,a)*Tx(2)+Nx(3,a)*Tx(3)) - N(a)*q)

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) 
     2         + wl*((N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b))
     3         + alp*(Nx(1,a)*Nx(1,b)+Nx(2,a)*Nx(2,b)+Nx(3,a)*Nx(3,b)))
         END DO
      END DO
      
      RETURN
      END SUBROUTINE eval3Adr
!---------------------------------------------------------------------
      PURE SUBROUTINE eval2Adr(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
      CLASS(adrType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN
      INTEGER, INTENT(IN) :: eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2   ks(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
      
      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,3D0,1D0/)
      INTEGER a, b, Ac
      REAL(KIND=8) T1, amd, wl, wr, Td, Tx(nsd), tauM, kU, 
     2   kssk, nTx, Tp, udTx, udNx(eNoN), u(nsd), C, alp, q

      T1  = eq%itg%af*eq%itg%gam*dt
      amd = eq%itg%am/T1
      C   = eq%mat%Cp()*eq%mat%rho()
      alp = eq%mat%K()/C
      q   = eq%mat%q()/C
      wr  = w*C
      wl  = wr*T1

      Td  = 0D0
      Tx  = 0D0
      DO a=1, eNoN
         Ac = eqN(a)
         Td = Td + N(a)*eq%T%A%s(Ac)
      
         Tx(1) = Tx(1) + Nx(1,a)*eq%T%s(Ac)
         Tx(2) = Tx(2) + Nx(2,a)*eq%T%s(Ac)
      END DO

      IF (ASSOCIATED(eq%U)) THEN
         u = 0D0
         DO a=1, eNoN
            Ac = eqN(a)
            u(1) = u(1) + N(a)*eq%U%v(1,Ac)
            u(2) = u(2) + N(a)*eq%U%v(2,Ac)
         END DO
 
         IF (ASSOCIATED(eq%dmn%Um)) THEN
            DO a=1, eNoN
               Ac = eqN(a)
               u(1) = u(1) - N(a)*eq%dmn%Um%v(1,Ac)
               u(2) = u(2) - N(a)*eq%dmn%Um%v(2,Ac)
            END DO
         END IF

         kU = u(1)*u(1)*ks(1,1) + u(2)*u(1)*ks(2,1) 
     2      + u(1)*u(2)*ks(1,2) + u(2)*u(2)*ks(2,2)

         kssk = ks(1,1)*ks(1,1) + ks(2,1)*ks(2,1) 
     2        + ks(1,2)*ks(1,2) + ks(2,2)*ks(2,2)
         
         nTx = ks(1,1)*Tx(1)*Tx(1) + ks(2,2)*Tx(2)*Tx(2)
     4       + (ks(1,2) + ks(2,1))*Tx(1)*Tx(2)
         IF (ISZERO(nTx)) nTx = eps

         udTx = u(1)*Tx(1)   + u(2)*Tx(2)
         udNx = u(1)*Nx(1,:) + u(2)*Nx(2,:)
         Tp   = ABS(Td + udTx)
         alp  = alp + Tp/SQRT(nTx)/2D0
         tauM = ct(4)/SQRT(ct(1)/dt/dt + ct(2)*kU + ct(3)*alp*alp*kssk)
         Tp   = -tauM*(Td + udTx)
      ELSE
         udTx = 0D0
         udNx = 0D0
         tauM = 0D0
      END IF

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + wr*(N(a)*(Td + udTx) - udNx(a)*Tp
     2      + alp*(Nx(1,a)*Tx(1)+Nx(2,a)*Tx(2)) - N(a)*q)
         
         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b)
     2         + wl*((N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)) 
     2         + alp*(Nx(1,a)*Nx(1,b)+Nx(2,a)*Nx(2,b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE eval2Adr
!---------------------------------------------------------------------
      PURE SUBROUTINE bEvalAdr(eq, eNoN, eqN, w, N, h, nV, lR, lK)
      CLASS(adrType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: eNoN, eqN(eNoN)
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN),
     2   lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)

      INTEGER a, b, Ac
      REAL(KIND=8) wl, T, udn

      wl = w*eq%itg%af*eq%itg%gam*dt
      
      IF (ASSOCIATED(eq%U)) THEN
         udn = 0D0
         T   = 0D0
         IF (ASSOCIATED(eq%dmn%Um)) THEN
            DO a=1, eNoN
               Ac = eqN(a)
               udn = udn + SUM(N(a)*(eq%U%v(:,Ac)-eq%dmn%Um%v(:,Ac))*nV)
               T   = T + N(a)*eq%T%s(Ac)
            END DO
         ELSE
            DO a=1, eNoN
               Ac  = eqN(a)
               udn = udn + SUM(N(a)*eq%U%v(:,Ac)*nV)
               T   = T + N(a)*eq%T%s(Ac)
            END DO
         END IF
         udn = 5D-1*(udn - ABS(udn))
      ELSE
         udn = 0D0
         T   = 0D0
      END IF

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*N(a)*(h(1) - udn*T)
         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) - wl*N(a)*N(b)*udn
         END DO
      END DO

      RETURN
      END SUBROUTINE bEvalAdr

!#####################################################################
!     To test the current implementation.
      SUBROUTINE TEST_ADRMOD(dmn, sct, gg, res)
      TYPE(dmnType), INTENT(IN) :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      TYPE(ggType), INTENT(IN) :: gg
      TYPE(varType), INTENT(INOUT) :: res

      INTEGER bType
      REAL(KIND=8) gs(1)
      TYPE(ggType) gl
      TYPE(adrType) eq

      eq = adrType(dmn, sct, 'Water', dmn%nMsh) ! ., ., mat_name, nBc
      IF (eq%mat%rho() .NE. 1D0) io%e = "Issue with newEq"
      io%o = "newEq: "//CLR("(PASSED)",3)
      bType = IBSET(0,bType_Dir)
      eq%bc(1) = bcType(dmn, "Face_1", gg, bType) ! ., fa_name, ., bType
      IF (dmn%nMsh .EQ. 2) THEN
         gs = (/1D-1/)
         gl = ggType(dmn,'Face_2',gs,0,0,0)
         CALL gl%cProfile(.FALSE.,.FALSE.,.TRUE.)
         CALL gl%map()
         eq%bc(2) = bcType(dmn, "Face_2", gl, bType)
      END IF
      io%o = "newBc: "//CLR("(PASSED)",3)
      CALL eq%ini()
      IF (ANY(eq%var(1)%A%sP .NE. 0D0)) io%e = "Issue with %setup"
      io%o = "eq%setup: "//CLR("(PASSED)",3)
      DO cTS=1, nTS
         CALL eq%update()
         IF (ANY(eq%T%n%s.NE.eq%T%o%s)) io%e="Issue with eq%update"
         IF (cTS .EQ. 1) io%o = "eq%update: "//CLR("(PASSED)",3)
         DO WHILE(.NOT.eq%satisfied())
            CALL eq%solve()
         END DO
         IF (cTS .EQ. 1) io%o = "eq%solve: "//CLR("(PASSED)",3)
         CALL eq%save()
         IF (cTS .EQ. 1) io%o = "eq%save: "//CLR("(PASSED)",3)
      END DO
      res = eq%T%varType

      RETURN
      END SUBROUTINE TEST_ADRMOD
      END MODULE ADRMOD
