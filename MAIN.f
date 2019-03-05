!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Main routine that contains the calls to all major routines and
!     general structure of the code.
!      
!--------------------------------------------------------------------

      PROGRAM MAIN
      USE AEQMOD
      USE VTKMOD
      IMPLICIT NONE

!     Major.minor version
      CHARACTER(LEN=*), PARAMETER :: version="1.0a"

!     MODULE VARIABLES
!     for communication between processors                  (CMMOD)
!     TYPE(cmType) :: cm 
!     Whether to use color in printing outputs              (IOMOD)
!     LOGICAL :: pClr = .TRUE.
!     Appended path to all files that are going to be saved (IOMOD)
!     CHARACTER(LEN=stdL) :: appPath = ""
!     A generic IO for my own use                           (IOMOD)
!     TYPE(ioType) :: io
!     Number of spatial dimensions                          (ELEMOD)
!     INTEGER :: nsd = 3
!     Number of time steps                                  (ITGMOD)
!     INTEGER nTS
!     Current time step                                     (ITGMOD)
!     INTEGER :: cTS = 0
!     Time step size                                        (ITGMOD)
!     REAL(KIND=8) dt
!     Physical time                                         (ITGMOD)
!     REAL(KIND=8) :: time = 0D0
!     Number of equations                                   (AEQMOD)
!     INTEGER nEq
!     All equations                                         (AEQMOD)
!     TYPE(cEqType), POINTER :: eq(:)

!     Following are the remaining variables with default value
!     LOGICAL VARIABLES
!     Solution strategy: coupled or uncoupled
      LOGICAL :: coupled = .TRUE.
!     Whether to averaged results
      LOGICAL :: saveAve = .FALSE.
!     Whether start from beginning or from a previous simulations
      LOGICAL :: iniFlag = .TRUE.
!     Whether to produce one restart file rather a series of them
      LOGICAL :: rstOW = .TRUE.
 
!     INTEGER VARIABLES
!     Number of initialization time steps
      INTEGER :: nITS = 0
!     Start saving after this number of time step
      INTEGER :: saveATS = 1
!     Increment in saving solutions
      INTEGER :: saveInc = 10
!     Increment in saving restart file
      INTEGER rstInc

!     CHARACTER VARIABLES      
!     Restart file name
      CHARACTER(LEN=stdL) :: rstName = "restart"
!     Stop_trigger file name
      CHARACTER(LEN=stdL) :: stopTrigName = "STOP_SIM"
!     Saved output file name
      CHARACTER(LEN=stdL) :: saveName = 'sol'

!     DERIVED TYPE VARIABLES
!     The structure to read fParam.mfs as the input file 
      TYPE(lstType) lst
!     The entire solution domain (meshes)
      TYPE(dmnType) :: dmn
!     To store the solution
      TYPE(vtkType) vtk
!     Initialization file
      TYPE(fileType) iniF
!     Restart file
      TYPE(fileType) rstF

!     And a few temporary variables
      LOGICAL ltmp, allOk, terminate
      INTEGER i, iEq, iMat, nMat
      REAL(KIND=8) idt
      CHARACTER(LEN=stdL) ctmp
      TYPE(matType) mat
      TYPE(lstType), POINTER :: lPtr
      TYPE(vtkType) vtmp

!!    Temporary Grant variable
      TYPE(prtType),POINTER :: prt
!---------------------------------------------------------------------
!     Starting the communicator
      cm = cmType()
      
!     Starting IO class to write to screen/file
      ctmp = "P#"//STR(cm%id())
      IF (cm%mas()) ctmp = ""
      io = ioTYpe(cm%mas(),tag=TRIM(ctmp))

!     Runing one of the test suites or read in mfs file
      i = IARGC()
      IF (i .EQ. 0) io%e = "Input .mfs file is missing as an argument"
      IF (i .GT. 1) io%e = "Too many arguments: only one is needed"
      CALL GETARG(1,ctmp)
      SELECT CASE (LOWER(ADJUSTL(ctmp)))
      CASE ('-ts1')
         CALL TEST_SUITE(1)
      CASE ('-ts2')
         CALL TEST_SUITE(2)
      CASE DEFAULT
         lst = lstType(ctmp)
      END SELECT

!     Adjusting root folder before anything else
      ltmp = .TRUE.
      lPtr => lst%get(ltmp,"Save results in a folder")
      IF (ltmp) appPath = STR(cm%np())//"-procs/"

!     Making fine adjustments to IO based on user request
      lPtr => lst%get(ltmp,"Verbose")
      IF (ASSOCIATED(lPtr)) CALL io%o%set(ltmp,ltmp)
      lPtr => lst%get(ltmp,"Warning")
      IF (ASSOCIATED(lPtr)) CALL io%w%set(ltmp,ltmp)
      lPtr => lst%get(ltmp,"Debug")
      IF (ASSOCIATED(lPtr)) CALL io%d%set(ltmp,ltmp)
      lPtr => lst%get(pClr,"Colorful terminal output")

!     Reading in all variables. Starting by the welcome message
      io%o = "                                              "
      io%o = " ----------------------------------------------"
      io%o = "             Copyright (C), 2013               "
      io%o = "                Mahdi Esmaily                  "
      io%o = "     Multi-Physics Finite Element Solver"
      io%o = CLR("            (MUPFES version "//version//")")
      io%o = " ----------------------------------------------"
      io%o = "                                               "
      io%o = " This program comes with ABSOLUTELY NO WARRANTY"
      io%o = " This is free software, and you are welcome to "
      io%o = " redistribute it under certain conditions."
      io%o = " For details see the LICENSE file."
      io%o = "                                               "
      io%o = " Number of processors: "//cm%np()
 
!     Solution control parameters 
      lPtr => lst%get(nsd,"Number of spatial dimensions",ll=2,ul=3)
      lPtr => lst%get(nITs,"Number of initialization time steps",ll=0)
      lPtr => lst%get(coupled,"Coupled solution strategy")
      lPtr => lst%get(nTs,"Number of time steps",1,ll=1)
      lPtr => lst%get(idt,"Time step size",1,lb=0D0)
      dt = idt

!     Reading the meshes through dmn
      dmn = dmnType(lst)
      CALL dmn%ini()

!     Triggering stop file 
      lPtr =>lst%get(stopTrigName,"Searched file name to trigger stop")
      stopTrigName = TRIM(appPath)//stopTrigName
      
!     Output vtk file
      lPtr => lst%get(saveName,"Name prefix of saved files")
      saveName = TRIM(appPath)//saveName
      vtk = vtkType(saveName,dmn)
      CALL vtk%open('w')
      CALL vtk%wPart()
      CALL vtk%close()
      lPtr => lst%get(saveInc,"Increment in saving files",ll=1)
      lPtr => lst%get(saveATS,"Start saving after time step",ll=1)
      lPtr => lst%get(saveAve,"Save averaged results")

!     Reading materials
      nMat = lst%srch("Add material",ll=1)
      DO iMat=1, nMat
         lPtr => lst%get(ctmp,"Add material",iMat)
         mat = matType(lPtr, ctmp)
         CALL ADD_MAT(mat)
      END DO

!     Reading all equations - implemented in AEQ 
      CALL READ_AEQ(dmn,lst)

!     Restart file
      rstInc  = saveInc
      lPtr   => lst%get(rstInc,"Increment in saving restart files",ll=1)
      lPtr   => lst%get(rstOW,"Overwrite restart file")
      lPtr   => lst%get(rstName,"Restart file name")
      rstName = TRIM(appPath)//rstName
      rstF    = fileType(TRIM(rstName)//"_last.bin",'binary')

!     Initializing simulation. Using restart file as default. 
      lPtr => lst%get(iniFlag,"Continue previous simulation")
      iniF = rstF 
      lPtr => lst%get(iniF,"Simulation initialization file path")
      IF (ASSOCIATED(lPtr) .AND. .NOT.iniFlag .AND. cm%mas()) io%w = 
     2   "Ignoring initialization file (simulation is restarted)"
      IF (iniFlag) THEN
         ctmp = iniF%name()
         i    = LEN(TRIM(ctmp))
         IF (ctmp(i-3:i) .EQ. ".bin") THEN
            IF (iniF%exist()) THEN
               CALL iniF%open('r')
               CALL RW_AEQ(dmn,iniF)
            END IF
         ELSE IF (ctmp(i-3:i) .EQ. ".vtk") THEN
            vtmp = vtkType(iniF%name(),dmn)
            IF (vtmp%exist()) THEN
               CALL vtmp%open('r')
               DO iEq=1, nEq
                  DO i=1, eq(iEq)%s%nVar
                     CALL vtmp%nRW(eq(iEq)%s%var(i)%n)
                  END DO
               END DO
               CALL vtmp%close()
            END IF
         ELSE
            io%e = "Unknown file extension for initialization"
         END IF
      END IF

!     Done with reading and initializing all variables.
      CALL lst%check()
      io%o = CLR(" Configuration is completed successfully",3)
      io%o = " "
      io%o = "---------------------------------------------------------"
      io%o = "Eq     N-i     T      dB   Ri/R0    R/Ri     lsIt  dB  %t"
      io%o = "---------------------------------------------------------"
!---------------------------------------------------------------------
!     Outer loop for marching in time.
      terminate = cTS .GE. nTS
      DO WHILE (.NOT.terminate)
!     Going to the next time step: updating variables
         DO iEq=1, nEq
            CALL eq(iEq)%s%update()
         END DO

         cTS = cTS + 1 ! Incrementing time step
         dt  = idt ! Calculating dt
         IF (cTS .LE. nITS) dt = idt/1D1 
         time = time + dt  ! Updating time

!     Solving all equation
         IF (coupled) THEN ! If equations are coupled
            DO ! non-linear iteration loop
               allOk = .TRUE.
               DO iEq=1, nEq
                  CALL eq(iEq)%s%solve()
                  allOk = allOk .AND. eq(iEq)%s%satisfied()
               END DO
               IF (allOk) EXIT
            END DO
         ELSE ! Uncoupled case
            DO iEq=1, nEq
               DO WHILE(.NOT.eq(iEq)%s%satisfied()) ! non-linear itr.
                  CALL eq(iEq)%s%solve()
               END DO
            END DO
         END IF

!     Writing csv flux/average files
         DO iEq=1, nEq
            CALL eq(iEq)%s%save()
         END DO

!     Deciding whether the simulation should be terminated
         IF (cm%mas()) INQUIRE(FILE=stopTrigName, EXIST=terminate)
         CALL cm%bcast(terminate) 
         terminate = terminate .OR. (cTS.GE.nTS)

!     Writing restart file
         IF (terminate .OR. MOD(cTS,rstInc).EQ.0) THEN
            IF (.NOT.rstOW) THEN
               CALL rstF%open('w',TRIM(rstName)//"_"//STR(cTS)//".bin")
            ELSE
               CALL rstF%open('w')
            END IF
            CALL RW_AEQ(dmn,rstF)
            IF (.NOT.rstOW .AND. cm%mas())
     2         CALL SYSTEM("ln -f "//TRIM(rstF%name())//" "//
     3         TRIM(rstName)//"_last.bin")
         END IF

!     Writing vtk file
         IF (cTS.GE.saveATS .AND. MOD(cTS,saveInc).EQ.0) THEN
            WRITE(ctmp,'(I3.3)') cTs
            IF (cTs .GE. 1000) ctmp = STR(cTs)
            CALL vtk%open('w',TRIM(saveName)//"_"//TRIM(ctmp))
            DO iEq=1, nEq
               DO i=1, eq(iEq)%s%nVar
                  CALL vtk%nRW(eq(iEq)%s%var(i)%n)
               END DO
            END DO
            CALL vtk%close()
         END IF
      END DO

!     Cleaning up (in reverse order of initialization) before exit 
      DO iEq=1, nEq
         CALL eq(iEq)%s%free()
      END DO
      DEALLOCATE(eq)
      CALL dmn%free()
      CALL cm%finalize()

      END PROGRAM MAIN
