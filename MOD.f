!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     All the data structures are defined in this module.
!      
!--------------------------------------------------------------------

      MODULE COMMOD
      USE CMMOD
      USE CHNLMOD
      

!--------------------------------------------------------------------      
!     Constants and upperbounds

!     maximum possible nsd, standard length for strings, 
!     random for history file handel, maximum possible number 
!     of output variables, size of blocks for openMP communications, 
!     master is assumed to have zero ID, version, maximum number of
!     properties, License expiration date
      INTEGER, PARAMETER :: maxnsd = 3, 

      CHARACTER, PARAMETER :: delimiter = "/"

!     Possible senarios for the output, followed by the possible outputs
!     Undefined output, extract it from A, extract it from Y, extract it
!     from D, calculate WSS, calculate vorticity, energy flux, heat
!     flux, absolute velocity (for FSI)
      INTEGER, PARAMETER :: outGrp_NA = 500, outGrp_A = 501,
     2   outGrp_Y = 502, outGrp_D = 503, outGrp_WSS = 504, 
     3   outGrp_vort = 505, outGrp_eFlx = 506, outGrp_hFlx = 507,
     4   outGrp_absV = 508, outGrp_DUInv = 509
      INTEGER, PARAMETER :: out_velocity = 599, out_pressure = 598,
     2   out_acceleration = 597, out_temperature = 596, out_WSS = 595,
     3   out_vorticity = 594, out_displacement = 593, 
     4   out_energyFlux = 592, out_heatFlux = 591, 
     5   out_absVelocity = 590, out_DUInv = 589


!--------------------------------------------------------------------
!     Here comes subTypes definitions later used in other derived types
!     Domain type is to keep track with element belong to which domain
!     and also different hysical quantities
      TYPE dmnType
!        The domain ID. Default includes entire domain
         INTEGER :: Id = -1
!        which physics must be solved in this domain
         INTEGER phys
!        The volume of this domain
         REAL(KIND=8) :: v = 0D0
!        physical properties, such as viscosity, density, ...
         REAL(KIND=8) :: prop(maxNProp) = 0D0
      END TYPE dmnType
     
!     Declared type for outputed variables
      TYPE outputType
!        Is this output suppose to be written into VTK, boundary, vol
         LOGICAL :: wtn(3) = .FALSE.
!        The group that this belong to (one of outType_*)
         INTEGER :: grp = outGrp_NA
!        Length of the outputed variable
         INTEGER l
!        Offset from the first index
         INTEGER o
!        The name to be used for the output and also in input file
         CHARACTER(LEN=stdL) name
      END TYPE outputType

     
!--------------------------------------------------------------------
!     All the subTypes are defined, now defining the major types that
!     will be directly allocated



      END MODULE COMMOD

!####################################################################
!     This subroutine is to read inputs of VTK files or boundaries.
!     nDOP(1) is the total number of outputs, nDOP(2) is the default 
!     number of outputs for VTK files, nDOP(3) is for boundaries, and
!     nDOP(4) is for volume
      SUBROUTINE READOUTPUTS(lEq, nDOP, outputs, list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER, INTENT(IN) :: nDOP(4), outputs(nDOP(1))
      TYPE(listType), INTENT(INOUT) :: list

      INTEGER nOut, iOut, i, j
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPO

      lEq%nOutput = nDOP(1)
      ALLOCATE(lEq%output(nDOP(1)))
      
      DO iOut=1, nDOP(1)
         SELECT CASE (outputs(iOut))
         CASE (out_velocity)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Velocity"
         CASE (out_pressure)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = nsd
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Pressure"
         CASE (out_acceleration)
            lEq%output(iOut)%grp  = outGrp_A
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Acceleration"
         CASE (out_temperature)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Temperature"
         CASE (out_displacement)
            lEq%output(iOut)%grp  = outGrp_D
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Displacement"
         CASE (out_WSS)
            lEq%output(iOut)%grp  = outGrp_WSS
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = maxnsd
            lEq%output(iOut)%name = "WSS"
         CASE (out_vorticity)
            lEq%output(iOut)%grp  = outGrp_vort
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = maxnsd
            lEq%output(iOut)%name = "Vorticity"
         CASE (out_energyFlux)
            lEq%output(iOut)%grp  = outGrp_eFlx
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Energy_flux"
         CASE (out_heatFlux)
            lEq%output(iOut)%grp  = outGrp_hFlx
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Heat_flux"
         CASE (out_absVelocity)
            lEq%output(iOut)%grp  = outGrp_absV
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Absolute_velocity"
         CASE (out_DUInv)
            lEq%output(iOut)%grp  = outGrp_DUInv
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "DU_invariants"
         CASE DEFAULT
            err = "Internal output undefined"
         END SELECT
      END DO

!     These are the default values, we use the first nDef/nBDef outputs
      DO j=1, 3
         lEq%output(1:nDOP(j+1))%wtn(j) = .TRUE.
      END DO
!     First reading the outputs for VTK files and then for boundaries 
!     and last for the volume
      nOut = list%srch("Output")
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output",iOut)
         SELECT CASE(TRIM(ctmp))
         CASE("Spatial")
            j = 1
         CASE("B_INT")
            j = 2
         CASE("V_INT")
            j = 3
         CASE DEFAULT
            j = -1
            err = TRIM(list%ping("Output",lPO))//" Undefined keyword"
         END SELECT
         DO i=1, lEq%nOutput
            lPtr => lPO%get(lEq%output(i)%wtn(j),lEq%output(i)%name)
         END DO
      END DO

      RETURN
      END SUBROUTINE READOUTPUTS

!####################################################################
!--------------------------------------------------------------------
!      
!     This routine is desinged to read the mesh/es that may come in
!     several formats and setup all parameters related to meshes
!      
!--------------------------------------------------------------------
!     The higher level routine that calls other routines      
      SUBROUTINE READMSH(list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list
 
      LOGICAL flag
      INTEGER i, j, iM, iFa, a, b, Ac, e
      REAL(KIND=8) maxX(nsd), minX(nsd), scaleF
      CHARACTER(LEN=stdL) ctmp, fName
      TYPE(listType), POINTER :: lPtr, lPM, lPBC
      TYPE(stackType) avNds
      TYPE(fileType) ftmp

      REAL(KIND=8), ALLOCATABLE :: tmpX(:,:), gX(:,:)

!     Setting dmnId parameter here, if there is at least one mesh that
!     has defined eId. 
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)
         
         lPtr => lPM%get(ftmp,"Domain file path")
         IF (ASSOCIATED(lPtr)) CALL msh(iM)%readEId(ftmp%open())
         
         lPtr => lPM%get(i,"Domain",ll=0,ul=BIT_SIZE(dmnId)-1)
         IF (ASSOCIATED(lPtr)) CALL msh(iM)%setEId(i)
         IF (ALLOCATED(msh(iM)%eId)) flag = .TRUE.
      END DO
      IF (flag) THEN
         ALLOCATE(dmnId(gtnNo))
         dmnId = 0
         DO iM=1, nMsh
            IF (.NOT.ALLOCATED(msh(iM)%eId)) CYCLE
            DO e=1, msh(iM)%gnEl
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%gIEN(a,e)
                  Ac = msh(iM)%gN(Ac)
                  dmnId(Ac) = IOR(dmnId(Ac),msh(iM)%eId(e))
               END DO
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE READMSH

!####################################################################
!     This function returns the domain that a set of nodes belongs to
      FUNCTION DOMAIN(lM, iEq, e)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iEq, e
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER DOMAIN

      INTEGER iDmn

      DOMAIN = 0
!     Domain Id of -1 counts for the entire domain      
      DO iDmn=1, eq(iEq)%nDmn
         DOMAIN = iDmn
         IF (eq(iEq)%dmn(iDmn)%Id .EQ. -1) RETURN
      END DO
      IF (.NOT. ALLOCATED(lM%eId)) err = "eId is not allocated,"//
     2   " while Id.NE.-1"
      DO iDmn=1, eq(iEq)%nDmn
         DOMAIN = iDmn
         IF (BTEST(lM%eId(e),eq(iEq)%dmn(iDmn)%Id)) RETURN
      END DO
      err = "Unable to find the domain ID of element "//e
      
      RETURN
      END FUNCTION DOMAIN

!--------------------------------------------------------------------
!     This function is true if "phys" is solved on "node" 
      PURE FUNCTION ISDOMAIN(iEq, node, phys)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iEq, node, phys
      LOGICAL ISDOMAIN

      INTEGER iDmn

      ISDOMAIN = .FALSE.
      IF (ALLOCATED(dmnId)) THEN
         DO iDmn=1, eq(iEq)%nDmn
            IF (eq(iEq)%dmn(iDmn)%phys .EQ. phys) THEN
               IF (BTEST(dmnId(node),eq(iEq)%dmn(iDmn)%Id)) THEN
                  ISDOMAIN = .TRUE.
                  RETURN
               END IF
            END IF
         END DO
      ELSE
!     No domain partitioning exists, so single domain is assumed and we
!     only need to check that
            IF (eq(iEq)%dmn(1)%phys .EQ. phys) ISDOMAIN = .TRUE.
      END IF

      RETURN
      END FUNCTION ISDOMAIN

!####################################################################

!     Distributing lM%dmnId if present to processors
      flag = ALLOCATED(dmnId)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(part(mg%gtnNo))
            part = dmnId
            DEALLOCATE(dmnId)
         ELSE
            ALLOCATE(part(0))
         END IF
         ALLOCATE(dmnId(mg%tnNo))
         dmnId = LOCAL(part)
         DEALLOCATE(part)
      END IF

!    Was a part of PICC - corrector
      IF (eq(cEq)%phys.EQ.phys_FSI) THEN
         s = eq(2)%s
         e = eq(2)%e
         DO Ac=1, tnNo
            IF (ISDOMAIN(cEq, Ac, phys_struct)) THEN
               An(s:e,Ac) = An(1:nsd,Ac)
               Yn(s:e,Ac) = Yn(1:nsd,Ac)
               Dn(s:e,Ac) = Dn(1:nsd,Ac)
            END IF
         END DO
      END IF

!        Assigns a value to eId array
         PROCEDURE :: assignEId => assignEIdMsh
!        Reads eId array from a file
         PROCEDURE :: readEId => readEIdMsh
!---------------------------------------------------------------------
!     Set domain ID to a given number for the entire or a range of
!     elements in a mesh      
      SUBROUTINE assignEIdMsh(lM, iDmn, ifirst, ilast)
      INTEGER, INTENT(IN) :: iDmn
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN), OPTIONAL :: ifirst, ilast
 
      INTEGER e, first, last

      first = 1
      IF (PRESENT(ifirst)) first = ifirst
      IF (first .LE. 0) io%e = "Out of range value in setDmnId"
      
      last = lM%nEl
      IF (PRESENT(ilast)) last = ilast
      IF (last.LT.first .OR. last.GT.lM%nEl) io%e = "Out of range"//
     2   " value in setDmnID"

      IF (.NOT.ALLOCATED(lM%eId)) THEN
         ALLOCATE(lM%eId(lM%nEl))
         lM%eId = 0
      END IF
      DO e=first, last
         lM%eId(e) = IBSET(lM%eId(e),iDmn)
      END DO

      RETURN
      END SUBROUTINE assignEIdMsh
!---------------------------------------------------------------------
!     Read domain from a file 
      SUBROUTINE readEIdMsh(lM, fid)
      CLASS(pMshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: fid
      
      INTEGER i, e, a
      CHARACTER(LEN=stdL) ctmp, stmp
      INTEGER, ALLOCATABLE :: iDmn(:)

      IF (.NOT.ALLOCATED(lM%eId)) THEN
         ALLOCATE(lM%eId(lM%nEl))
         lM%eId = 0
      END IF

!     Check to see if I need to increase the size of eId while
!     changing it
      ALLOCATE (iDmn(BIT_SIZE(lM%eId)))
      ctmp = "Domain ID exceeds BIT_SIZE(eId). Reduce the number"//
     2   " of domains or increase BIT_SIZE(eId)."
      DO e=1, lM%nEl
         READ (fid,"(A)") stmp
         i = CheckNoNumbers(stmp)
         IF (i .GT. BIT_SIZE(lM%eId)) io%e = ctmp
         READ (stmp,*) iDmn(1:i)
         DO a=1, i
            IF (iDmn(a) .GE. BIT_SIZE(lM%eId)) THEN
               io%e = ctmp
            ELSE IF (iDmn(a) .LT. 0) THEN
               io%e = "Domain ID must greater or equal to 0"
            END IF
            lM%eId(e) = IBSET(lM%eId(e),iDmn(a))
         END DO
      END DO
      CLOSE (fid)

      RETURN
      END SUBROUTINE readEIdMsh

      IF (ALLOCATED(this%eId)) DEALLOCATE(this%eId)

!     This it to distribute eId, if allocated
         flag = .FALSE.
         IF (ALLOCATED(lM%eId)) THEN
            flag = .TRUE.
            ALLOCATE(part(lM%gnEl))
            DO e=1, lM%gnEl
               Ec = otnPtr(e)
               part(Ec) = lM%eId(e)
            END DO
            DEALLOCATE(lM%eId)
         END IF

      CALL cm%bcast(flag)

!     Communicating eId, if neccessary
      IF (flag) THEN
         ALLOCATE(lM%eId(nEl))
         CALL cm%scatter(part,lM%eId)
         IF (ALLOCATED(part)) DEALLOCATE(part)
      END IF

!!!!!!!!! from BAFAINI
      IF (mvMsh) THEN
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct)) i = i + 1
         END DO
         ALLOCATE(gNodes(i))
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct)) THEN
               i = i + 1
               gNodes(i) = a
            END IF
         END DO
         CALL memLS_BC_CREATE(lhs, nFacesLS, i, nsd, BC_TYPE_Dir,gNodes)
      END IF

