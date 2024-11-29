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
FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
FFLAGS += -J$(MOD_DIR)
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
	@echo "  done cleaning"

cleanall : clean
	rm -fr OBJ/obj*
	rm -f lib*.a
	@echo "  done all cleaning"
#===============================================
#============= Main executable and tests  ======
#===============================================
test_cavity_mode.exe: $(OBJ_DIR)/test_cavity_mode.o $(LIBA)
	$(FFC) -o test_cavity_mode.exe  $(FFLAGS) $(OBJ_DIR)/test_cavity_mode.o $(LIBA)

$(OBJ_DIR)/test_cavity_mode.o:$(TESTS_DIR)/test_cavity_mode.f90
	$(FFC) -c -o $(OBJ_DIR)/test_cavity_mode.o $(FFLAGS) $(TESTS_DIR)/test_cavity_mode.f90
#=================================================================================
#=================================================================================
$(OBJ_DIR)/MC_cavity_mode_m.o:$(SRC_DIR)/MC_cavity_mode_m.f90
	$(FFC) -c -o $(OBJ_DIR)/MC_cavity_mode_m.o $(FFLAGS) $(SRC_DIR)/MC_cavity_mode_m.f90
$(OBJ_DIR)/MC_operator_1D_m.o:$(SRC_DIR)/MC_operator_1D_m.f90
	$(FFC) -c -o $(OBJ_DIR)/MC_operator_1D_m.o $(FFLAGS) $(SRC_DIR)/MC_operator_1D_m.f90
#=================================================================================
#============= module dependencies =============
#===============================================
$(OBJ_DIR)/test_cavity_mode.o:       $(LIBA) 


