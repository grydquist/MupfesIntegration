!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     Contains an abstract class for Equation. Implements calls that are
!     common among all equations: solve, stat, const, and etc. Other
!     equations must provide %eval implementation.
!     For usage and testing eq class, see TEST_???MOD with ??? being
!     equations that use eq class. 
!      
!---------------------------------------------------------------------
      MODULE EQMOD
      USE GGMOD
      USE SCTMOD
      USE CBCMOD
      USE LSMOD
      USE MATMOD
      IMPLICIT NONE

!     boundary conditions types. Items of eq list can be combined:
!     Drichlet, Neumann, periodic, coupled, resistance, backflow
!     stabilization, impose BC on the integral of state variable or D
!     (instead of Y), diplacement dependent
      INTEGER, PARAMETER :: bType_Dir = 0, bType_Neu = 1, bType_per = 2,
     3   bType_cpl = 3, bType_res = 4, bType_bfs = 5, bType_impD = 6, 
     4   bType_ddep = 7

!     Possible formulation for equation: based on spatial configuration,
!     based on displacement at last time step, based on reference
!     configuration
      INTEGER, PARAMETER :: eqConf_NA = 200, eqConf_spt = 201, 
     2   eqConf_old = 202, eqConf_ref = 203

!     Volumetric condition data type
      TYPE vcType
         PRIVATE
!        Pointer to memLS%bc
         INTEGER lsPtr
!        On which variable eq BC should be imposed on
         INTEGER :: iVar = 1
!        Whether VC is imposed. will have dimension [dmn%dof]
         LOGICAL, ALLOCATABLE :: impose(:)
!        The variable pointing to the imposed values
         TYPE(gVarType), POINTER :: g => NULL()
      END TYPE vcType

      INTERFACE vcType
         PROCEDURE :: newVc
      END INTERFACE vcType
!---------------------------------------------------------------------
!     Boundary condition data type
      TYPE bcType
         PRIVATE
!        The face_ID associated with this BC
         INTEGER iFa 
!        The mesh_ID associated with this BC
         INTEGER iM
!        Pre/Res/Flat/Para... boundary types
         INTEGER :: bType = 0
!        Direction for imposing the BC
         INTEGER :: eDrn = 0
!        Pointer to memLS%bc
         INTEGER lsPtr
!        On which variable eq BC should be imposed on
         INTEGER :: iVar = 1
!        Neu: defined resistance
         REAL(KIND=8) :: r = 0D0
!        The face corresponding to this BC
         TYPE(faceType), POINTER :: fa => NULL()
!        The mesh corresponding to this BC
         TYPE(mshType), POINTER :: lM => NULL()
!        Prescribed values at bundaries (For either Dir or Neu)
         TYPE(ggType) :: gg
      CONTAINS 
!        Deallocates all variables
         PROCEDURE, PRIVATE :: free => freeBc
      END TYPE bcType

      INTERFACE bcType
         PROCEDURE :: newBc, newBcFl
      END INTERFACE bcType
!---------------------------------------------------------------------
!     A container for equations
      TYPE :: cEqType
         CLASS(eqType), POINTER :: s => NULL()
      END TYPE cEqType
!---------------------------------------------------------------------
!     Equation type. Extending solution control 
      TYPE, EXTENDS(sctType), ABSTRACT :: eqType
!        Whether the equation is initialized
         LOGICAL :: izd = .FALSE.
!        Number of BCs
         INTEGER :: nBc = 0
!        Number of tracked variables
         INTEGER :: nVar = 0
!        Number of sub equations
         INTEGER :: nSub = 0
!        Configuration for calculation of shape functions
         INTEGER :: conf = eqConf_spt
!        Equation symbol (to be overwritten)
         CHARACTER(LEN=2) :: sym = "NA"
!        Time integrator
         TYPE(itgType) :: itg
!        Linear solver
         TYPE(lsType) :: ls
!        For coupled BCs
         TYPE(cbcType) :: cbc
!        Properties are stored in this variable
         TYPE(matType), POINTER :: mat => NULL()
!        Pointer to the domain where eq equation is solved
         TYPE(dmnType), POINTER :: dmn => NULL()
!        Which sub equation to solve for a given mesh
         INTEGER, ALLOCATABLE :: subPtr(:)
!        VCs associated with this equation
         TYPE(vcType) :: vc
!        BCs associated with this equation
         TYPE(bcType), ALLOCATABLE :: bc(:)
!        All variables tracked by eq equation
         TYPE(gVarType), POINTER :: var(:)
!        Files for saving boundary avergaes/fluxes
         TYPE(fileType), ALLOCATABLE :: bAF(:)
!        If this equation has sub equations (e.g. FSI)
         TYPE(cEqType), ALLOCATABLE :: sub(:)
      CONTAINS
!        Solves the governing equation
         PROCEDURE :: solve => solveEq
!        Prints stats of eq equation into screen/txt-files
         PROCEDURE :: stat => statEq
!        Updating solution, old=new
         PROCEDURE :: update => updateEq
!        Initializing an equation: to be called once
         PROCEDURE :: ini => iniEq
!        Write the boundary average/fluxes in txt files
         PROCEDURE :: save => saveEq
!        Deallocates all variables
         PROCEDURE :: free => freeEq
!        Constructs a part of the structure
         PROCEDURE :: newEq
         PROCEDURE :: newEqFl
         GENERIC :: new => newEq, newEqFl
!        Constructs the local stiffness matrix and assembels it
         PROCEDURE, PRIVATE :: cnst => cnstEq
!        Sets Dirichlet BC
         PROCEDURE, PRIVATE :: setDirBc => setDirBcEq
!        Sets VC
         PROCEDURE, PRIVATE :: setVc => setVcEq
!        Initializing a boundary condition
         PROCEDURE, PRIVATE :: iniBc => iniBcEq
!        Initializing a volumetric condition
         PROCEDURE, PRIVATE :: iniVc => iniVcEq
!        Constructs the local stiffness matrix and assembels it
         PROCEDURE, PRIVATE :: cnstBc => cnstBcEq
!        Evaulate the governing equation in 3D and 2D
         PROCEDURE(evalEqIf), DEFERRED :: eval3, eval2
!        Evaulate the governing equation boundary contribution
         PROCEDURE(bEvalEqIf), DEFERRED :: bEval
!        Setting up equation-specific parameters
         PROCEDURE(setupEqIf), DEFERRED :: setup
      END TYPE eqType

      INTERFACE
         SUBROUTINE setupEqIf(eq, var)
            IMPORT
            CLASS(eqType), INTENT(INOUT), TARGET :: eq
            TYPE(gVarType), INTENT(IN), TARGET, OPTIONAL :: var(:)
         END SUBROUTINE setupEqIf

         SUBROUTINE evalEqIf(eq, eNoN, eqN, w, J, N, Nx, ks, lR, lK)
            IMPORT
            CLASS(eqType), INTENT(IN) :: eq
            INTEGER, INTENT(IN) :: eNoN
            INTEGER, INTENT(IN) :: eqN(eNoN)
            REAL(KIND=8), INTENT(IN) :: w, J, N(eNoN), Nx(nsd,eNoN),
     2         ks(nsd,nsd)
            REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2         lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
         END SUBROUTINE evalEqIf
         
         SUBROUTINE bEvalEqIf(eq, eNoN, eqN, w, N, h, nV, lR, lK)
            IMPORT
            CLASS(eqType), INTENT(IN) :: eq
            INTEGER, INTENT(IN) :: eNoN
            INTEGER, INTENT(IN) :: eqN(eNoN)
            REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h(:), nV(nsd)
            REAL(KIND=8), INTENT(INOUT) :: lR(eq%ls%nU,eNoN), 
     2         lK(eq%ls%nU*eq%ls%nU,eNoN,eNoN)
         END SUBROUTINE 
      END INTERFACE 
!---------------------------------------------------------------------
!     Parameters that are equations-specific, to be provided when
!     creating a new eq
      TYPE eqSpType
         INTEGER nVar
         INTEGER itgO
         INTEGER lsAl
         CHARACTER(LEN=2) sym
      END TYPE eqSpType
!---------------------------------------------------------------------
!     To keep track of time since start of simulation
      REAL(KIND=8), PRIVATE :: time_keeper(3) = 0D0

      CONTAINS
!#####################################################################

!     IMPLEMENTATION

!#####################################################################
      FUNCTION newBc(dmn, name, gg, bType, eDrn, r, cbc) RESULT(bc)
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(ggType), INTENT(IN) :: gg
      INTEGER, INTENT(IN) :: bType
      INTEGER, INTENT(IN), OPTIONAL :: eDrn
      REAL(KIND=8), INTENT(IN), OPTIONAL :: r
      TYPE(cbcType), INTENT(INOUT), OPTIONAL :: cbc
      TYPE(bcType) :: bc
      
!     Finding msh that eq BC belongs to
      CALL dmn%find(name, bc%iM, bc%iFa)
      bc%fa   => dmn%msh(bc%iM)%fa(bc%iFa)
      bc%lM   => dmn%msh(bc%iM)
      bc%gg    = gg
      bc%bType = bType
      IF (PRESENT(eDrn)) bc%eDrn  = eDrn
      IF (PRESENT(r))    bc%r     = r
      IF (BTEST(bType,bType_cpl)) THEN
         IF (.NOT.PRESENT(cbc)) io%e = "newBc: missing cbc argument"
         IF (BTEST(bType,bType_Dir)) CALL cbc%addBc(name, cbc_Dir)
         IF (BTEST(bType,bType_Neu)) CALL cbc%addBc(name, cbc_Neu)
      END IF
      IF (ANY(dmn%prjFa.EQ.name) .AND. cm%mas()) io%w = "Make sure "//
     1   "one/both faces are prescribed as Neu/Dir, respectively, "//
     2   "when imposing BC on a projected face."

      RETURN
      END FUNCTION newBc
!---------------------------------------------------------------------
      FUNCTION newBcFl(lst, dmn, name, cbc) RESULT(bc)
      TYPE(lstType), INTENT(INOUT) :: lst
      TYPE(dmnType), INTENT(IN) :: dmn
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(cbcType), INTENT(INOUT), OPTIONAL :: cbc
      TYPE(bcType) :: bc
      
      LOGICAL ltmp
      INTEGER bType, eDrn
      REAL(KIND=8) r
      CHARACTER(LEN=stdL) ctmp
      TYPE(lstType), POINTER :: lPtr
      TYPE(ggType) gg

!     Default values
      bType = 0
      eDrn  = 0
      r     = 0D0

!     Reading the type: Dir/Neu/Per
      lPtr => lst%get(ctmp,"Type")
      SELECT CASE (ctmp)
      CASE ("dirichlet","dir")
         bType = IBSET(bType,bType_Dir)
      CASE ("neumann","neu")
         bType = IBSET(bType,bType_Neu)
      CASE ("periodic","per")
         bType = IBSET(bType,bType_per)
         io%e = "Periodic BC hasn't been implemented yet"
      CASE DEFAULT
         io%e = TRIM(lst%ping("Type",lPtr))//" Unexpected BC type"
      END SELECT

!     Impose BC on the state variable or its integral
      ltmp = .FALSE.
      lPtr => lst%get(ltmp,"Impose on state variable integral")
      IF (ltmp) bType = IBSET(bType,bType_impD)

!     Limit BC to an specific direction
      lPtr => lst%get(eDrn,"Effective direction")

      ctmp = "Steady"
      lPtr => lst%get(ctmp,"Time dependency")
      SELECT CASE (ctmp)
      CASE ('coupled')
         bType = IBSET(bType,bType_cpl)
      CASE ('resistance')
         bType = IBSET(bType,bType_res)
         IF (.NOT.BTEST(bType,bType_Neu)) io%e = "Resistance"//
     2      " is only defined for Neu BC"
         lPtr => lst%get(r,"Value",1)
      END SELECT

!     Constructing imposed profile
      gg = ggType(lst, dmn, name)
      IF (BTEST(bType,bType_Dir)) CALL gg%map()

!     Finally constructing the BC
      bc = newBc(dmn, name, gg, bType, eDrn, r, cbc)

      RETURN
      END FUNCTION newBcFl
!--------------------------------------------------------------------
      FUNCTION newVc(g, impose) RESULT(vc)
      TYPE(gVarType), INTENT(IN), TARGET :: g
      LOGICAL, INTENT(IN) :: impose(g%dmn%nMsh)
      TYPE(vcType) :: vc
      
      vc%g => g
      ALLOCATE(vc%impose, SOURCE=impose)

      RETURN
      END FUNCTION newVc
!---------------------------------------------------------------------
      SUBROUTINE iniBcEq(eq, iBc, lsPtr)
      CLASS(eqType), INTENT(INOUT) :: eq
      INTEGER, INTENT(IN) :: iBc
      INTEGER, INTENT(INOUT) :: lsPtr

      INTEGER a, e, Ac, g
      REAL(KIND=8) n(nsd)
      INTEGER, ALLOCATABLE :: gN(:)
      REAL(KIND=8), ALLOCATABLE :: sV(:,:), sVl(:,:)
      ASSOCIATE(bc => eq%bc(iBc), lM => eq%bc(iBc)%lM, dmn => eq%dmn, 
     2   fa => eq%bc(iBc)%fa, iFa => eq%bc(iBc)%iFa)
       
      IF (BTEST(bc%bType,bType_Dir)) THEN
         ALLOCATE(sVl(nsd,fa%dof), sV(nsd,lM%nNo))
         lsPtr    = lsPtr + 1
         bc%lsPtr = lsPtr
         IF (bc%eDrn .NE. 0) THEN
            sVl = 1D0
            sVl(bc%eDrn,:) = 0D0
         ELSE
            sVl = 0D0
         END IF
         CALL eq%ls%addBc(lsPtr, fa%dof, BC_TYPE_Dir, fa%uP, sVl)
      ELSE IF (BTEST(bc%bType,bType_Neu)) THEN
         ALLOCATE(sVl(nsd,fa%nNo), sV(nsd,lM%nNo), gN(fa%nNo))
         DO a=1, fa%nNo
            gN(a) = lM%uP(fa%gN(a))
         END DO
         IF (BTEST(bc%bType,bType_res).OR.BTEST(bc%bType,bType_cpl))THEN
            sV = 0D0
            DO e=1, fa%nEl
               DO g=1, fa%nG
                  n = lM%normal(iFa,e,g)
                  DO a=1, fa%eNoN
                     Ac = fa%IEN(a,e)
                     sV(:,Ac) = sV(:,Ac) + fa%N(a,g)*fa%w(g)*n
                  END DO
               END DO
            END DO
            DO a=1, fa%nNo
               Ac       = fa%gN(a)
               sVl(:,a) = sV(:,Ac)
            END DO
            lsPtr    = lsPtr + 1
            bc%lsPtr = lsPtr
            CALL eq%ls%addBc(lsPtr, fa%nNo, BC_TYPE_Neu, gN, sVl)
         ELSE
            bc%lsPtr = 0
         END IF
      ELSE
         io%e = "iniBc: Unxpected bType"
      END IF

      RETURN
      END ASSOCIATE
      END SUBROUTINE iniBcEq
!--------------------------------------------------------------------
      SUBROUTINE iniVcEq(eq, lsPtr)
      CLASS(eqType), INTENT(INOUT) :: eq
      INTEGER, INTENT(INOUT) :: lsPtr

      INTEGER a, Ac, nNo, iM, i, dof
      INTEGER, ALLOCATABLE :: impose(:), gImpose(:), gN(:)
      REAL(KIND=8), ALLOCATABLE :: sVl(:,:)

      dof = eq%vc%g%dof
      IF (eq%var(eq%vc%iVar)%dof .NE. dof) io%e = 
     2   "iniVcEq: incompatible dof between var and g"

!     Finding the nodes that this volumetric condition is imposed on
      ALLOCATE(impose(eq%dmn%dof), gImpose(eq%dmn%gdof))
      impose = 0
      DO iM=1, eq%dmn%nMsh
         IF (.NOT.eq%vc%impose(iM)) CYCLE
         DO a=1, eq%dmn%msh(iM)%nNo
            i = eq%dmn%msh(iM)%uP(a)
            impose(i) = 1
         END DO
      END DO

      gImpose = 0
      WHERE(impose.EQ.1)
         gImpose(eq%dmn%guP) = 1
      END WHERE
      gImpose = cm%reduce(gImpose)
      DEALLOCATE(eq%vc%impose)
      ALLOCATE(eq%vc%impose(eq%dmn%dof))
      
      eq%vc%impose = .FALSE.
      nNo = COUNT(gImpose.GT.0)
      ALLOCATE(sVl(dof,nNo), gN(nNo))
      i = 0
      DO a=1, eq%dmn%dof
         Ac = eq%dmn%guP(a)
         IF (gImpose(Ac) .EQ. 1) THEN
            i = i + 1
            gN(i) = a
            eq%vc%impose(a) = .TRUE.
         END IF
      END DO
      sVl         = 0D0
      lsPtr       = lsPtr + 1
      eq%vc%lsPtr = lsPtr
      CALL eq%ls%addBc(lsPtr, nNo, BC_TYPE_Dir, gN, sVl)
      DEALLOCATE(impose, gImpose, gN)

      RETURN
      END SUBROUTINE iniVcEq
!---------------------------------------------------------------------
      SUBROUTINE setDirBcEq(eq, iBc)
      CLASS(eqType), INTENT(IN) :: eq
      INTEGER, INTENT(IN) :: iBc
      
      INTEGER a, Ac, i, lb, ub, lI, pI
      REAL(KIND=8), POINTER :: pA(:), pY(:)
      REAL(KIND=8), ALLOCATABLE :: lA(:,:), lY(:,:)
      ASSOCIATE(bc => eq%bc(iBc), fa => eq%bc(iBc)%fa, 
     2   var => eq%var(eq%bc(iBc)%iVar), gg => eq%bc(iBc)%gg)

      IF (BTEST(bc%bType,bType_impD)) THEN
         pA => var%n%sP
         pY => var%Dn%sP
      ELSE
         pA => var%An%sP
         pY => var%n%sP
      END IF

      IF (bc%eDrn .EQ. 0) THEN
         lb = 1
         ub = var%dof
      ELSE
         lb = bc%eDrn
         ub = bc%eDrn
      END IF
      
      i = ub - lb + 1
      ALLOCATE(lA(i,fa%dof), lY(i,fa%dof))
      CALL gg%eval(time, lY, lA)

      DO a=1, fa%dof
         Ac = fa%uP(a)
         DO i=lb, ub
            lI = i - lb + 1
            pI = (Ac-1)*var%dof + i
            pY(pI) = lY(lI,a)
            pA(pI) = lA(lI,a)
         END DO
      END DO

      RETURN
      END ASSOCIATE
      END SUBROUTINE setDirBcEq
!---------------------------------------------------------------------
      SUBROUTINE setVcEq(eq)
      CLASS(eqType), INTENT(IN) :: eq
      
      INTEGER a
      ASSOCIATE(var => eq%var(eq%vc%iVar))

      IF (var%dof .EQ. 1) THEN
         WHERE (eq%vc%impose)
            var%Dn%s = eq%vc%g%Dn%s
            var% n%s = eq%vc%g% n%s
            var%An%s = eq%vc%g%An%s
         END WHERE
      ELSE IF (var%dof .EQ. nsd) THEN
         DO a=1, eq%dmn%dof
            IF (eq%vc%impose(a)) THEN
               var%Dn%v(:,a) = eq%vc%g%Dn%v(:,a)
               var% n%v(:,a) = eq%vc%g% n%v(:,a)
               var%An%v(:,a) = eq%vc%g%An%v(:,a)
            END IF
         END DO
      ELSE
         io%e = "setVcEq: undefined dof"
      END IF

      RETURN
      END ASSOCIATE
      END SUBROUTINE setVcEq
!---------------------------------------------------------------------      
!     eq routine calculates the parameters that are required for
!     sub-equations and calculates them based on general variables.
      SUBROUTINE cnstBcEq(eq, iBc)
      CLASS(eqType), INTENT(INOUT), TARGET :: eq
      INTEGER, INTENT(IN) :: iBc
      
      INTEGER a, g, eNoN, nU, e, Ac, i
      REAL(KIND=8) w, nV(nsd), Jac, Q
      REAL(KIND=8), ALLOCATABLE :: lK(:,:,:), lR(:,:), N(:), lV(:), 
     2   h(:), hl(:,:), hg(:,:)
      INTEGER, ALLOCATABLE :: eqN(:)
      CLASS(eqType), POINTER :: lEq
      ASSOCIATE(bc => eq%bc(iBc), iFa => bc%iFa, fa => bc%fa,
     2   lM => bc%lM, var => eq%var(bc%iVar), dof => var%dof)

!     Figuring out which equation to solve
      lEq => eq
      IF (eq%nSub .NE. 0) lEq => eq%sub(eq%subPtr(bc%iM))%s

      eNoN = fa%eNoN
      nU   = lEq%ls%nU
      ALLOCATE(lK(nU*nU,eNoN,eNoN), lR(nU,eNoN), N(eNoN), eqN(eNoN), 
     2   lV(lM%nNo), h(dof), hl(dof,eNoN), hg(dof,lM%nNo))
      
!     Since eDrn might be nonzero
      hg = 0D0 
      CALL bc%gg%eval(time, hg)
      IF (BTEST(bc%bType,bType_res)) THEN
         IF (var%dof .EQ. 1) THEN
            DO a=1, fa%nNo
               Ac = fa%gN(a)
               i  = lM%uP(Ac)
               lV(Ac) = var%s(i)
            END DO
         ELSE IF (var%dof .EQ. nsd) THEN
            DO a=1, fa%nNo
               Ac = fa%gN(a)
               i  = lM%uP(Ac)
               lV(Ac) = DOT_PRODUCT(var%v(:,i),fa%nV(:,a))
            END DO
         ELSE
            io%e = "cnstBc: unknown var%dof for bType_res"
         END IF
         Q  = lM%integ(iFa,lV)
         Q  = cm%reduce(Q)
         hg = (Q*bc%r)*hg
      END IF
 
      IF (BTEST(bc%bType,bType_ddep)) THEN
         IF (var%dof .EQ. 1) THEN
            DO a=1, fa%nNo
               Ac = fa%gN(a)
               i  = lM%uP(Ac)
               hg(:,Ac) = hg(:,Ac)*var%D%s(i)
            END DO
         ELSE IF (var%dof .EQ. nsd) THEN
            DO a=1, fa%nNo
               Ac = fa%gN(a)
               i  = lM%uP(Ac)
               hg(:,Ac) = hg(:,Ac)*DOT_PRODUCT(var%D%v(:,i),fa%nV(:,a))
            END DO
         ELSE
            io%e = "cnstBc: unknown var%dof for bType_ddep"
         END IF
      END IF

      DO e=1, fa%nEl
         DO a=1, eNoN
            Ac      = fa%IEN(a,e)
            eqN(a)  = lM%uP(Ac)
            hl(:,a) = hg(:,Ac)
         END DO
         lK = 0D0
         lR = 0D0
         DO g=1, fa%nG
            nV  = lM%normal(iFa, e, g)
            Jac = NORM2(nV)
            nV  = nV/Jac
            w   = fa%w(g)*Jac
            N   = fa%N(:,g)
            h   = 0D0
            DO a=1, eNoN
               h = h + N(a)*hl(:,a)
            END DO
            
            CALL lEq%bEval(eNoN, eqN, w, N, h, nV, lR, lK)
         END DO
         CALL eq%ls%assem(eNoN, eqN, lK, lR)
      END DO
      DEALLOCATE(eqN, lK, lR, N, lV)

      RETURN
      END ASSOCIATE
      END SUBROUTINE cnstBcEq
!---------------------------------------------------------------------
      SUBROUTINE freeBc(bc)
      CLASS(bcType), INTENT(INOUT) :: bc

      CALL bc%gg%free()

      bc%bType = 0
      bc%eDrn  = 0
      bc%iVar  = 1
      bc%r     = 0D0
      bc%fa   => NULL()
      bc%lM   => NULL()

      RETURN 
      END SUBROUTINE freeBc
!#####################################################################
      SUBROUTINE newEq(eq, eqSp, dmn, sct, mName, nBc, ls, itg, cbc)
      CLASS(eqType), INTENT(INOUT), TARGET :: eq
      TYPE(eqSpType), INTENT(IN) :: eqSp
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      TYPE(sctType), INTENT(IN) :: sct
      CHARACTER(LEN=*), INTENT(IN) :: mName
      INTEGER, INTENT(IN) :: nBc
      TYPE(lsType), INTENT(IN), OPTIONAL :: ls
      TYPE(itgType), INTENT(IN), OPTIONAL :: itg
      TYPE(cbcType), INTENT(IN), OPTIONAL :: cbc

      eq%nVar    = eqSp%nVar
      eq%sym     = eqSp%sym
      eq%sctType = sct
      eq%nBc     = nBc
      eq%mat    => FIND_MAT(mName)
      eq%dmn    => dmn

      IF (PRESENT(ls)) THEN
         IF (ls%isIzd()) io%e = "newEq: ls should not be initialized"
         eq%ls = ls
      ELSE
         eq%ls = lsType(eqSp%lsAl, nBc+1, dmn) ! +1 for vc
      END IF
      IF (PRESENT(itg)) THEN
         eq%itg = itg
      ELSE
         eq%itg = itgType(eqSp%itgO)
      END IF
      IF (PRESENT(cbc)) eq%cbc = cbc
      ALLOCATE(eq%bc(eq%nBc), eq%var(eqSp%nVar), eq%bAF(eqSp%nVar))
      IF (eq%nVar .LE. 0) io%e = eq%sym//": nVar <= 0"
      IF (eq%nBc  .LT. 0) io%e = eq%sym//": nBc < 0"
      IF (.NOT.ASSOCIATED(eq%mat)) io%e = eq%sym//
     2   "Material was not found"
      IF (.NOT.ASSOCIATED(eq%dmn)) io%e = eq%sym//
     2   "Domain is undefined"
      
      RETURN
      END SUBROUTINE newEq     
!---------------------------------------------------------------------
      SUBROUTINE newEqFl(eq, eqSp, dmn, lst)
      CLASS(eqType), INTENT(INOUT), TARGET :: eq
      TYPE(eqSpType), INTENT(IN) :: eqSp
      TYPE(dmnType), INTENT(IN), TARGET :: dmn
      TYPE(lstType), INTENT(INOUT) :: lst

      INTEGER nBc, iBc
      CHARACTER(LEN=stdL) ctmp
      TYPE(lstType), POINTER :: lPtr
      TYPE(sctType) sct
      TYPE(lsType) ls
      TYPE(cbcType) cbc

!     Reading the type of material and number of BCs
      lPtr => lst%get(ctmp,"Material",1)
      nBc = lst%srch("Add BC")
!     Reading sct, ls, and cbc from lst and calling the default constructor
      sct = sctType(lst)
      ls  = lsType(lst, eqSp%lsAl, nBc+1, dmn)
      cbc = cbcType(lst)
      CALL newEq(eq, eqSp, dmn, sct, ctmp, nBc, ls=ls, cbc=cbc)

!     Adjusting symbol if neccessary
      lPtr => lst%get(ctmp,"Symbol")
      IF (ASSOCIATED(lPtr)) eq%sym = ctmp(1:LEN(eq%sym))

!     Searching for BCs
      io%o = " Number of imposed BC for equation <"//TRIM(eq%sym)//">: "
     2   //eq%nBc
      DO iBc=1, eq%nBc
         lPtr => lst%get(ctmp,"Add BC",iBc)
         eq%bc(iBc) = bcType(lPtr, dmn, ctmp, eq%cbc)
      END DO

      RETURN
      END SUBROUTINE newEqFl 
!---------------------------------------------------------------------
      SUBROUTINE iniEq(eq)
      CLASS(eqType), INTENT(INOUT), TARGET :: eq
      
      INTEGER iBc, lsPtr, iVar, nU
      CHARACTER(LEN=stdL) ctmp

      IF (.NOT.eq%dmn%izd) io%e = "iniEq: Un-initialized domain"
      IF (eq%izd) THEN
         io%w = "Equation <"//eq%sym//"> already initialzied"
         RETURN
      END IF
      eq%izd = .TRUE.

!     Equation-specific setup of the variables
      CALL eq%setup()

!     Initializing linear solver with nU DOF per node.
      nU = SUM(eq%var%dof)
      CALL eq%ls%ini(nU)

      lsPtr = 0
!     Initializing BCs and imposing Dirichlet BCs
      DO iBc=1, eq%nBc
         IF (.NOT.ASSOCIATED(eq%bc(iBc)%fa)) 
     2      io%e = "iniEq: BC <"//iBc//"> is not assigned"
         CALL eq%iniBc(iBc,lsPtr)
         IF (BTEST(eq%bc(iBc)%bType,bType_Dir)) CALL eq%setDirBc(iBc)
      END DO
!     Initializing VCs and imposing them
      IF (ASSOCIATED(eq%vc%g)) THEN
         CALL eq%iniVc(lsPtr)
         CALL eq%setVc()
      END IF

!     Prepapring the boundary files 
      DO iVar=1, eq%nVar
         ctmp = TRIM(appPath)//eq%sym//"_"//TRIM(eq%var(iVar)%name)//
     2      "_integ.csv"
         eq%bAF(iVar) = fileType(ctmp,'binary')
         CALL eq%bAF(iVar)%open('a')
         IF (cm%slv()) CALL eq%bAF(iVar)%close()
      END DO

!     Creating csv headers
      CALL eq%save()

!     Starting the timer
      time_keeper(1) = CPUT()

      RETURN
      END SUBROUTINE iniEq
!---------------------------------------------------------------------
!     Add contribution of current equation to the LHS/RHS         
      SUBROUTINE cnstEq(eq, iM)
      CLASS(eqType), INTENT(INOUT), TARGET :: eq
      INTEGER, INTENT(IN) :: iM
      
      INTEGER g, eNoN, e, nU, a, Ac
      REAL(KIND=8) w, J, ks(nsd,nsd)
      INTEGER, ALLOCATABLE :: eqN(:)
      REAL(KIND=8), ALLOCATABLE :: lK(:,:,:), lR(:,:), N(:), Nx(:,:),
     2   d(:,:), xl(:,:)
      CLASS(eqType), POINTER :: lEq
      ASSOCIATE(dmn => eq%dmn, lM => dmn%msh(iM))

      IF (.NOT.eq%izd) io%e = "constEq: Un-initialized equation"

!     Figuring out which equation to solve
      lEq => eq
      IF (eq%nSub .NE. 0) lEq => eq%sub(eq%subPtr(iM))%s

      eNoN = lM%eNoN
      nU   = lEq%ls%nU
      ALLOCATE(d(nsd,eNoN), xl(nsd,eNoN), eqN(eNoN), 
     2   lK(nU*nU,eNoN,eNoN), lR(nU,eNoN), Nx(nsd,eNoN), N(eNoN))

      DO e=1, lM%nEl
!     Setting intial values
         lK = 0D0
         lR = 0D0

!     Change reference configuration when mesh is moving
         IF (ASSOCIATED(dmn%Um)) THEN
            IF (lEq%conf .EQ. eqConf_spt) THEN
               DO a=1, eNoN
                  Ac      = lM%IEN(a,e)
                  eqN(a)  = lM%uP(Ac)
                  xl(:,a) = lM%x(:,Ac) + dmn%Um%D %v(:,lM%uP(Ac))
               END DO
            ELSE IF (lEq%conf .EQ. eqConf_old) THEN
               DO a=1, eNoN
                  Ac      = lM%IEN(a,e)
                  eqN(a)  = lM%uP(Ac)
                  xl(:,a) = lM%x(:,Ac) + dmn%Um%Do%v(:,lM%uP(Ac))
               END DO
            ELSE IF (lEq%conf .EQ. eqConf_ref) THEN
               DO a=1, eNoN
                  Ac      = lM%IEN(a,e)
                  eqN(a)  = lM%uP(Ac)
                  xl(:,a) = lM%x(:,Ac)
               END DO
            ELSE
               io%e = "cnstEq: Undefined eq%conf"
            END IF
         ELSE
            DO a=1, eNoN
               Ac      = lM%IEN(a,e)
               eqN(a)  = lM%uP(Ac)
               xl(:,a) = lM%x(:,Ac)
            END DO
         END IF

!     Computing shape functions and looping over Gauss's points
         DO g=1, lM%nG
            IF (g.EQ.1.OR..NOT.lM%lShpF) CALL lM%dNdx(g, xl, Nx, J, ks)
            IF (ISZERO(J)) io%e = "Jac < 0 @ element "//e
            w = lM%w(g)*J
            N = lM%N(:,g)
            IF (nsd .EQ. 3) THEN
               CALL lEq%eval3(eNoN, eqN, w, J, N, Nx, ks, lR, lK)
            ELSE IF (nsd .EQ. 2) THEN
               CALL lEq%eval2(eNoN, eqN, w, J, N, Nx, ks, lR, lK)
            ELSE
               io%e = "solveEq: undefined %eval for this nsd"
            END IF
         END DO

!     Now doing the assembly part
         CALL eq%ls%assem(eNoN, eqN, lK, lR)
      END DO
      DEALLOCATE(d, xl, eqN, lK, lR, Nx, N)
      io%d = "Mesh <"//TRIM(lM%name)//"> is assembled"

      RETURN
      END ASSOCIATE
      END SUBROUTINE cnstEq
!---------------------------------------------------------------------
      SUBROUTINE solveEq(eq)
      CLASS(eqType), INTENT(INOUT) :: eq

      INTEGER iM, iVar, iBc, i
      INTEGER, ALLOCATABLE :: incL(:)
      REAL(KIND=8), ALLOCATABLE :: res(:)

      IF (.NOT.eq%izd) io%e = "solveEq: Un-initialized equation"

      ALLOCATE(res(eq%ls%nFa), incL(eq%ls%nFa))
!     Solving 0D domain equations, if needed
      IF (eq%cbc%isIzd()) THEN
         CALL eq%cbc%solve(eq%var(1), eq%var(2))
         DO iBc=1, eq%nBc
            IF (BTEST(eq%bc(iBc)%bType,bType_cpl)) CALL eq%cbc%get(
     2         eq%bc(iBc)%fa%name, eq%bc(iBc)%gg%gs, eq%bc(iBc)%r)
         END DO
      END IF
      
!     Setting R and Val to zero
      CALL eq%ls%reset()

!     Imposing Dirichlet BCs
      DO iBc=1, eq%nBc
         IF (BTEST(eq%bc(iBc)%bType,bType_Dir)) CALL eq%setDirBc(iBc)
      END DO
      IF (ASSOCIATED(eq%vc%g)) CALL eq%setVc()

!     Initiator step
      DO iVar=1, eq%nVar
         CALL eq%itg%initiate(eq%var(iVar))
      END DO

      io%d = "Assembling equation in volume <"//eq%sym//">"
      DO iM=1, eq%dmn%nMsh
         CALL eq%cnst(iM)
      END DO

!     Constructing the element stiffness matrix for boundaries
      incL = 0
      DO iBc=1, eq%nBc
         IF (BTEST(eq%bc(iBc)%bType,bType_Neu)) CALL eq%cnstBc(iBc)
         i = eq%bc(iBc)%lsPtr
         IF (i .NE. 0) THEN
            res(i)  = eq%itg%gam*dt*eq%bc(iBc)%r
            incL(i) = 1
         END IF
      END DO
      IF (ASSOCIATED(eq%vc%g)) incL(eq%vc%lsPtr) = 1

      io%d = "Solving linear system for equation <"//eq%sym//">"
      CALL eq%ls%solve(incL, res)

!     Updating norm for solution control
      CALL eq%upNorm(eq%ls%RI%iNorm)

!     Solution is obtained, now updating (Corrector)
      CALL eq%itg%correct(eq%var, eq%ls%soln())
      
!     Checking for exceptions
      CALL io%w%checkException()

!     Writing out the time passed, residual, and etc.
      CALL eq%stat()

      RETURN
      END SUBROUTINE solveEq
!---------------------------------------------------------------------
!     Updates the solution
      SUBROUTINE updateEq(eq)
      CLASS(eqType), INTENT(INOUT) :: eq
      
      INTEGER iVar

      DO iVar=1, eq%nVar
         eq%var(iVar)%Ao%sP = eq%var(iVar)%An%sP
         eq%var(iVar)% o%sP = eq%var(iVar)% n%sP
         eq%var(iVar)%Do%sP = eq%var(iVar)%Dn%sP
         CALL eq%itg%predict(eq%var(iVar))
      END DO
      eq%itr = 0
 
      RETURN  
      END SUBROUTINE updateEq
!---------------------------------------------------------------------
!     Prepares the output of mupfes to the standard output.
      SUBROUTINE statEq(eq)
      CLASS(eqType), INTENT(INOUT) :: eq

      INTEGER i
      REAL(KIND=8) tmp, tmp2
      CHARACTER c1, c2
      CHARACTER(LEN=stdL) sOut

      IF (cm%slv()) RETURN
      
      time_keeper(3) = CPUT() - time_keeper(1)
      sOut = eq%sym//" "//STR(cTS,5)//"-"//ADJUSTL(STR(eq%itr,2))//
     2   " "//STR(time_keeper(3),6)

      IF (ISZERO(eq%zNorm)) THEN
         tmp  = 1D0
         tmp2 = 1D0
         i    = 0
      ELSE
         tmp  = eq%ls%RI%iNorm/MAX(eq%zNorm,eps)
         tmp2 = eq%ls%RI%fNorm/MAX(eq%ls%RI%iNorm,eps)
         i    = INT(2D1*LOG10(tmp/MAX(eq%pNorm,eps)))
      END IF

      IF (i .GT. 20) THEN
         c1 = "!"; c2 = "!"
      ELSE
         c1 = "["; c2 = "]"
      END IF
      sOut = TRIM(sOut)//"  "//c1//STR(i,4)//" "//STR(tmp,7)//" "//
     2   STR(tmp2,7)//c2

      IF (ISZERO(time_keeper(3),time_keeper(2))) time_keeper(3) =
     2   (1D0+eps)*time_keeper(2) + eps

      tmp = 1D2*eq%ls%RI%callD/(time_keeper(3) - time_keeper(2))
      time_keeper(2) = time_keeper(3)
      IF (ABS(tmp) .GT. 1D2) tmp = 1D2

      IF (eq%ls%RI%suc) THEN
         c1 = "["; c2 = "]"
      ELSE
         c1 = "!"; c2 = "!"
      END IF
      io%o = TRIM(sOut)//"  "//c1//STR(eq%ls%RI%itr,4)//" "//
     2   STR(NINT(eq%ls%RI%dB),3)//" "//STR(NINT(tmp),3)//c2

      RETURN
      END SUBROUTINE statEq
!---------------------------------------------------------------------
      SUBROUTINE saveEq(eq)
      CLASS(eqType), INTENT(INOUT) :: eq

      INTEGER, PARAMETER :: pcn = 13
      INTEGER iM, iFa, iVar, pos
      CHARACTER(LEN=(SUM(eq%dmn%msh%nFa)+eq%dmn%nMsh+1)*(pcn+2)-1) s, t
      REAL(KIND=8), ALLOCATABLE :: rS(:), rV(:,:)
      TYPE(mshType), POINTER :: lM

      IF (.NOT.eq%izd) io%e = "saveEq: Un-initialized equation"

!     This is the special case that is called only once
      IF (cTS .EQ. 0) THEN
         s   = "mesh/face"
         t   = "volume/area  , "//STR(eq%dmn%msh%vol,pcn)
         pos = pcn + 1
         DO iM=1, eq%dmn%nMsh
            t = TRIM(t)//", "//STR(eq%dmn%msh(iM)%fa%area,pcn)
            s(pos:pos+pcn+1) = ", "//eq%dmn%msh(iM)%name(1:pcn)
            pos = pos + pcn + 2
         END DO
         DO iM=1, eq%dmn%nMsh
            DO iFa=1, eq%dmn%msh(iM)%nFa
               s(pos:pos+pcn+1)=", "//eq%dmn%msh(iM)%fa(iFa)%name(1:pcn)
               pos = pos + pcn + 2
            END DO
         END DO
         s(pos:pos) = eol
         t(pos:pos) = eol
         DO iVar=1, eq%nVar
            CALL eq%bAF(iVar)%setPos(1)
            CALL eq%bAF(iVar)%rw(s)
            CALL eq%bAF(iVar)%rw(t)
         END DO
         RETURN
      END IF
!     In general, however:
      DO iVar=1, eq%nVar
         s   = STR(cTS)
         pos = pcn + 1
         DO iM=1, eq%dmn%nMsh
            lM => eq%dmn%msh(iM)
            IF (eq%var(iVar)%dof .EQ. 1) THEN
               ALLOCATE(rS(lM%nNo))
               rS = eq%var(iVar)%n%s(lM%uP)
               s(pos:pos+pcn+1) = ", "//STR(lM%integ(rS),pcn)
               pos = pos + pcn + 2
               DO iFa=1, lM%nFa
                  s(pos:pos+pcn+1) = ", "//STR(lM%integ(iFa,rS),pcn)
                  pos = pos + pcn + 2
               END DO
               DEALLOCATE(rS)
            ELSE IF (eq%var(iVar)%dof .EQ. nsd) THEN
               ALLOCATE(rV(nsd,lM%nNo))
               rV = eq%var(iVar)%n%v(:,lM%uP)
               s(pos:pos+pcn+1) = ", "//STR(lM%integ(rV),pcn)
               pos = pos + pcn + 2
               DO iFa=1, lM%nFa
                  s(pos:pos+pcn+1) = ", "//STR(lM%integ(iFa,rV),pcn)
                  pos = pos + pcn + 2
               END DO
               DEALLOCATE(rV)
            ELSE
               io%e = "saveEq: undefined var%dof"
            END IF
         END DO
         s(pos:pos) = eol
         CALL eq%bAF(iVar)%setPos(LEN(s)*(1+cTS)+1)
         CALL eq%bAF(iVar)%rw(s)
      END DO

      RETURN 
      END SUBROUTINE saveEq
!---------------------------------------------------------------------
      SUBROUTINE freeEq(eq)
      CLASS(eqType), INTENT(INOUT) :: eq

      INTEGER iBc, iVar

      CALL eq%ls%free()

      IF (ALLOCATED(eq%bc)) THEN
         DO iBc=1, eq%nBc
            CALL eq%bc(iBc)%free()
         END DO
         DEALLOCATE(eq%bc)
      END IF
      IF (ASSOCIATED(eq%var)) THEN
         DO iVar=1, eq%nVar
            CALL eq%bAF(iVar)%close()
            CALL eq%var(iVar)%free()
         END DO
         DEALLOCATE(eq%var, eq%bAF)
      END IF
      eq%izd  = .FALSE.
      eq%nBc  = 0
      eq%nVar = 0
      eq%dmn => NULL()
      eq%mat => NULL()

      RETURN 
      END SUBROUTINE freeEq

      END MODULE EQMOD
