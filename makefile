#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
 FFC = gfortran

 
 
#=================================================================================
OBJ_DIR=OBJ/obj
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
#
MOD_DIR=$(OBJ_DIR)
EXTMod= 

SRC_DIR=SRC
MAIN_DIR=APP
TESTS_DIR=TESTS

LIB=libMolecCav
LIBA=$(LIB).a
#=================================================================================
# External Libraries directory
OOPT=0
OOMP=1
LLAPACK=0
INT=4
ext_obj=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)

ExtLibDIR := Ext_Lib
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))
#
QD_DIR=$(ExtLibDIR)/QDUtilLib
QDMOD_DIR=$(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA=$(QD_DIR)/libQD$(ext_obj).a
#===============================================================================
EXTLib     = $(QDLIBA)
EXTMod     = -I$(QDMOD_DIR)

#=================================================================================
FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
FFLAGS += -J$(MOD_DIR) $(EXTMod)
#=================================================================================
#=================================================================================
SRCFILES=MC_cavity_mode_m.f90 MC_operator_1D_m.f90

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))

$(info ***********************************************************************)
$(info ***********************************************************************)
$(info ***********COMPILER:     $(FFC))
$(info ***********OBJ:          $(OBJ))
$(info ***********FFLAGS:       $(FFLAGS))
$(info ***********************************************************************)
$(info ***********************************************************************)

#===============================================
#============= Tests ===========================
#===============================================
.PHONY: ut UT
UT ut: test_cavity_mode.exe
	./test_cavity_mode.exe < data.nml > test.log
	grep "successfully" test.log
	@echo "  done Tests"

#===============================================
#============= all: lib, tests ...  ============
#===============================================
.PHONY: all
all: $(LIBA) test_cavity_mode.exe
#===============================================
#============= Library: libQD.a  ===============
#===============================================
.PHONY: lib
lib: $(LIBA)

$(LIBA): $(OBJ)
	ar -cr $(LIBA) $(OBJ)
	@echo "  done Library: "$(LIBA)
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f test*.exe
	rm -f test.log
	@echo "  done cleaning"

cleanall : clean
	rm -fr OBJ/obj*
	rm -f lib*.a
	cd Ext_Lib ; ./cleanlib
	@echo "  done all cleaning"
#===============================================
#============= Main executable and tests  ======
#===============================================
test_cavity_mode.exe: $(OBJ_DIR)/test_cavity_mode.o $(LIBA)
	$(FFC) -o test_cavity_mode.exe  $(FFLAGS) $(OBJ_DIR)/test_cavity_mode.o $(LIBA) $(EXTLib)

$(OBJ_DIR)/test_cavity_mode.o:$(TESTS_DIR)/test_cavity_mode.f90
	$(FFC) -c -o $(OBJ_DIR)/test_cavity_mode.o $(FFLAGS) $(TESTS_DIR)/test_cavity_mode.f90
#=================================================================================
#=================================================================================
$(OBJ_DIR)/MC_cavity_mode_m.o:$(SRC_DIR)/MC_cavity_mode_m.f90
	$(FFC) -c -o $(OBJ_DIR)/MC_cavity_mode_m.o $(FFLAGS) $(SRC_DIR)/MC_cavity_mode_m.f90
$(OBJ_DIR)/MC_operator_1D_m.o:$(SRC_DIR)/MC_operator_1D_m.f90
	$(FFC) -c -o $(OBJ_DIR)/MC_operator_1D_m.o $(FFLAGS) $(SRC_DIR)/MC_operator_1D_m.f90
#=================================================================================
#===============================================
#===============================================
#== external libraries
#
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib
$(QDLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA)
#
#===============================================
#============= module dependencies =============
#===============================================
$(OBJ_DIR)/test_cavity_mode.o:       $(LIBA)

$(OBJ):            | $(QDLIBA)


