#     Copyright, 2013
#     Mahdi Esmaily Moghadam

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#--------------------------------------------------------------------
#
#     This is the Makefile to build mupfes.
#
#--------------------------------------------------------------------

# HERE COMES THE DEFINITIONS
# PRT is below INS below
TOP = /mnt/f/GradSchool/mupfes/trunk
include $(TOP)/Makefile.in

BIN_DIR      = $(TOP)/bin
PARMETIS_DIR = $(TOP)/external/ParMetis-3.2.0
memLS_DIR    = $(TOP)/memLS
cplBC_DIR    = $(TOP)/cplBC
MPI_DIR      = $(TOP)/external/dummy_mpi

memLS_INC = -I$(memLS_DIR)/include
memLS_LIB = -L$(memLS_DIR) -lmemLS
METIS_INC = -I$(PARMETIS_DIR)/ParMETISLib
METIS_LIB = -L$(PARMETIS_DIR) -lparmetis -lmetis -lm
cplBC_INC = -I$(cplBC_DIR)/include
MPI_INC   = -I$(MPI_DIR)
MPI_LIB   = -L$(MPI_DIR) -lmpi
FFLAGS	  = -g -check all

MYEXEC = $(BIN_DIR)/mupfes.1.0a
MYFUN  = UTIL.f \
         COMU.f \
         IO.f \
         FILE.f \
         LST.f \
         MAT.f \
         ELE.f \
         PMSH.f \
         MSH.f \
         VAR.f \
         VTK.f \
         GG.f \
         ITG.f \
         SCT.f \
         CBC.f \
         LS.f \
         EQ.f \
         ADR.f \
         LED.f \
         INS.f \
		 PRT.f \
         SVK.f \
         FSI.f \
         BBO.f \
         AEQ.f \
         TEST.f \
         MAIN.f \

OTHER	=  POST.f \
         MOD.f \
         LOADNRB.f \
         NURBS.f \
         PRT.f \

ifeq ($(seq),1)
   INCLUDES = $(memLS_INC) $(MPI_INC) $(cplBC_INC) 
   LIBS     = $(memLS_LIB) $(MPI_LIB)
else
   INCLUDES = $(memLS_INC) $(cplBC_INC)
   LIBS     = $(METIS_LIB) $(memLS_LIB)
   MYFUN   += SPLIT.f
endif

#############################################################################
# AND HERE ARE THE RULES

SRC = $(patsubst %,$(SRC_DIR)/%,$(MYFUN)) 
OBJ = $(patsubst %.f,$(OBJ_DIR)/%.o,$(MYFUN)) 

.PHONY: $(MYEXEC)
$(MYEXEC): $(OBJ) $(memLS_DIR)/libmemLS.a $(BIN_DIR)
	$(FORTRAN) $(OBJ) $(LIBS) $(FFLAGS) -lc -o $@

$(OBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/UTIL.o : $(SRC_DIR)/UTIL.f  
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/COMU.o : $(SRC_DIR)/COMU.f
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/IO.o : $(SRC_DIR)/IO.f $(OBJ_DIR)/UTIL.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/FILE.o : $(SRC_DIR)/FILE.f $(OBJ_DIR)/UTIL.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/LST.o : $(SRC_DIR)/LST.f $(OBJ_DIR)/IO.o $(OBJ_DIR)/FILE.o $(OBJ_DIR)/COMU.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/MAT.o : $(SRC_DIR)/MAT.f $(OBJ_DIR)/LST.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/ELE.o : $(SRC_DIR)/ELE.f $(OBJ_DIR)/IO.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/PMSH.o : $(SRC_DIR)/PMSH.f $(OBJ_DIR)/ELE.o $(OBJ_DIR)/LST.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/MSH.o : $(SRC_DIR)/MSH.f $(OBJ_DIR)/PMSH.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/VAR.o : $(SRC_DIR)/VAR.f $(OBJ_DIR)/MSH.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/VTK.o : $(SRC_DIR)/VTK.f $(OBJ_DIR)/VAR.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/GG.o : $(SRC_DIR)/GG.f $(OBJ_DIR)/VAR.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/ITG.o : $(SRC_DIR)/ITG.f $(OBJ_DIR)/VAR.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/SCT.o : $(SRC_DIR)/SCT.f $(OBJ_DIR)/LST.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/CBC.o : $(SRC_DIR)/CBC.f $(OBJ_DIR)/ITG.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/LS.o : $(SRC_DIR)/LS.f $(OBJ_DIR)/VAR.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/EQ.o : $(SRC_DIR)/EQ.f $(OBJ_DIR)/GG.o $(OBJ_DIR)/SCT.o $(OBJ_DIR)/CBC.o $(OBJ_DIR)/LS.o $(OBJ_DIR)/MAT.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/ADR.o : $(SRC_DIR)/ADR.f $(OBJ_DIR)/EQ.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/LED.o : $(SRC_DIR)/LED.f $(OBJ_DIR)/EQ.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/INS.o : $(SRC_DIR)/INS.f $(OBJ_DIR)/EQ.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@
	
# Here's PRT.f
$(OBJ_DIR)/PRT.o : $(SRC_DIR)/PRT.f $(OBJ_DIR)/EQ.o   
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/SVK.o : $(SRC_DIR)/SVK.f $(OBJ_DIR)/EQ.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/FSI.o : $(SRC_DIR)/FSI.f $(OBJ_DIR)/INS.o $(OBJ_DIR)/SVK.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/BBO.o : $(SRC_DIR)/BBO.f $(OBJ_DIR)/INS.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/AEQ.o : $(SRC_DIR)/AEQ.f $(OBJ_DIR)/ADR.o $(OBJ_DIR)/LED.o $(OBJ_DIR)/FSI.o $(OBJ_DIR)/BBO.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/TEST.o : $(SRC_DIR)/TEST.f $(OBJ_DIR)/VTK.o $(OBJ_DIR)/AEQ.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/MAIN.o : $(SRC_DIR)/MAIN.f $(OBJ_DIR)/VTK.o $(OBJ_DIR)/AEQ.o
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) $(METIS_INC) -c $< -o $@

clean:
	rm -r -f $(OBJ_DIR) $(MYEXEC) *.mod

