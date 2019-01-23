!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!---------------------------------------------------------------------
!      
!     To run all the tests of all classes. 
!      
!---------------------------------------------------------------------
      SUBROUTINE TEST_SUITE(nMsh)
      USE AEQMOD
      USE VTKMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nMsh

      REAL(KIND=8), PARAMETER :: 
     1   adr_bm(2) = (/0.275875038856942D0,0.245447045878021D0/),
     2   led_bm(2) = (/-0.66075899086489D0,-0.36911618980853D0/),
     3   ins_bm(2) = (/-0.07688622630753D0,-0.06522031151993D0/), 
     4   cbc_bm(2) = (/-0.01268848955696D0,-6.4883285314866D-3/),
     5   svk_bm(2) = (/-0.75162801886320D0,-0.41887974096508D0/),
     5   fsi_bm(2) = (/-0.75162801886320D0,-0.01558553520764D0/),
     5   bbo_bm(2) = (/-0.07718366326774D0,-6.5375395538401D-2/)
      TYPE(dmnType) dmn
      TYPE(varType) prs, vel
      TYPE(gVarType) gPrs
      TYPE(ggType) gg
      TYPE(sctType) sct
      TYPE(insType) ns
      TYPE(varType) res

!     Following are module variables
      cm  = cmType()
      io  = ioType(cm%mas())
      nsd = 3
      dt  = 1D-1
      nTS = 2

      CALL HDR('COMU')
      CALL TEST_CMMOD()
      CALL HDR('IO')
      CALL TEST_IOMOD()
      CALL HDR('FILE')
      IF (cm%mas()) CALL TEST_FILEMOD()
      CALL HDR('LST')
      CALL TEST_LSTMOD()
      CALL HDR('MAT')
      CALL TEST_MATMOD()
      CALL HDR('ELE')
      CALL TEST_ELEMOD()
      CALL HDR('MSH')
      CALL TEST_MSHMOD()
      CALL HDR('VAR')
      CALL TEST_VARMOD(nMsh, dmn, prs, vel, gPrs)
      CALL HDR('VTK')
      CALL TEST_VTKMOD(prs, vel)
      CALL HDR('GG')
      CALL TEST_GGMOD(dmn, gg)
      CALL HDR('ITG')
      CALL TEST_ITGMOD(gPrs)
      CALL HDR('SCT')
      CALL TEST_SCTMOD(sct)
      CALL HDR('LS')
      CALL TEST_LSMOD(dmn)
      res = varType(1,'dummy_var',dmn)
      CALL HDR('ADR')
      CALL TEST_ADRMOD(dmn, sct, gg, res)
      CALL VRFY(res,adr_bm)
      CALL HDR('LED')
      CALL TEST_LEDMOD(dmn, sct, gg, res)
      CALL VRFY(res,led_bm)
      CALL HDR('INS')
      CALL TEST_INSMOD(dmn, sct, gg, ns, res)
      CALL VRFY(res,ins_bm)
      CALL HDR('SVK')
      CALL TEST_SVKMOD(dmn, sct, gg, res)
      CALL VRFY(res,svk_bm)
      CALL HDR('FSI')
      CALL TEST_FSIMOD(dmn, sct, gg, res)
      CALL VRFY(res,fsi_bm)
      CALL HDR('BBO')
      CALL TEST_BBOMOD(ns, sct, gg, res)
      CALL VRFY(res,bbo_bm)
      CALL HDR('CBC')
      CALL TEST_CBCMOD(dmn, sct, gg, res)
      CALL VRFY(res,cbc_bm)
       
      CALL gg%free()
      CALL prs%free()
      CALL vel%free()
      CALL gPrs%free()
      CALL dmn%free()
      CALL cm%finalize()

      STOP
      CONTAINS 
!---------------------------------------------------------------------
      SUBROUTINE HDR(s)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: s
      
      CALL cm%barrier()
      io%o = ""
      io%o = "-----------------------------------------------"
      io%o = "          < < Testing "//s//" module > > "
      io%o = "-----------------------------------------------"

      RETURN
      END SUBROUTINE HDR
!---------------------------------------------------------------------
      SUBROUTINE VRFY(u, eu)
      IMPLICIT NONE
      TYPE(varType), INTENT(IN) :: u
      REAL(KIND=8), INTENT(IN) :: eu(2)

      INTEGER, PARAMETER :: Ac = 1
      REAL(KIND=8), PARAMETER :: tol = 1D-6
      REAL(KIND=8), ALLOCATABLE :: Ug(:), Ul(:)
      ASSOCIATE(lM => u%dmn%msh(1))
      
      ALLOCATE(Ul(lM%nNo), Ug(lM%gnNo))
      Ul = u%s(lM%uP)
      CALL cm%global(Ul, Ug, lM%gN)
      IF (cm%mas() .AND. NINT(Ug(Ac)/eu(dmn%nMsh)/tol).NE.NINT(1D0/tol))
     2   io%e = "Result invalid: "//STR(Ug(Ac),19)//" .NE. "//
     3   STR(eu(dmn%nMsh),19)
      io%o = "Verifying against benchmark: "//CLR("(PASSED)",3)
      DEALLOCATE(Ug,Ul)

      RETURN
      END ASSOCIATE
      END SUBROUTINE VRFY
      END SUBROUTINE TEST_SUITE
