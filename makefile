##################################################################################
##################################################################################
##################################################################################
# execute this file by command line using the following syntaxes :
# "make" to execute the first command defined in this file below
# "make <command>" to execute a particular command defined below
# "make -f <name of this file> <command (optional)>" if another file in Make is... 
# ...present in this directory
##################################################################################
##################################################################################
##################################################################################

#=================================================================================
#========================Definition of the Frotran Compiler=======================
#=================================================================================
# /!\ /!\ /!\. 
# in a variable definition, the value of the var MUST not be followed by ANY space (otherwise it is counted in the string and brongs havoc) 
# /!\ /!\ /!\.
FFC = gfortran
# syntaxe to create a variable (FFC) in Make, and attribute to the string "gfortran". NB: FC = Fortran Compiler
# the name of the fortran compiler to use here  

# Possible values: (Empty: gfortran)
#                  gfortran (version: 9.0 linux and osx)


#=================================================================================
#==================Definition of other useful names in variables==================
#=============================1. Of the main librairy=============================
#=================================================================================
MAIN_DIR  = APP
# the name of the directory where to find the source file(s) of the application/exemple program
MAIN      = App_MolecCav
# the name of the application file without the extension
OBJ_DIR   = OBJ/obj
# the name of the directory where to store the objects .o and .mod files (/obj because we might... 
# ...want to create different sub libraries with different parameters => different obj/ directories)
MOD_DIR   = $(OBJ_DIR)
# in make the variables values are treated and accessed the same way as in shell/bash : between "$(<name_of_the_variable>)"
# the name of the directory where to store the .mod files (here same as OBJ but can be different)
SRC_DIR   = SRC
# the name of the directory where to find the source files of the library's modules
TESTS_DIR = TESTS
# the name of the directory where to find the source files of the library's test programs
OUTPUT_DIR = OUT
# the name of the directory where to store the log files of the library's programs
DATA_DIR = DATA
# the name of the directory where to find the data files of the library

EXTMod    =                      
# the directory where the external libraries are to find (cf next subpart)

LIB      = libMolecCav
# the name of the library (the actual name of the actual library, not the name of any file of the library so far)
LIBA     = $(LIB).a
# the name (/!\ WITH the extension) of the library (the name of its full static object file this time = all... 
# ...the .o files of its modules)

$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)) # the instruction A = $(shell <shell command>) in Make will store into the var A the...
                                                  # ...result of the execution of the <shell command> command in the shell (and thereby...
                                                  # ...execute this command in the shell). Here the result of the command is note stored...
                                                  # ...but still executed.
                                                  # as we did in the script, this more compact command test the existence of OBJ/ and creates...
                                                  # ...it if it was not here (-p allow to create also the sub-directory /obj/)

$(shell [ -d $(OUTPUT_DIR) ] || mkdir -p $(OUTPUT_DIR))


#=================================================================================
#==================Definition of other useful names in variables==================
#==========================2. Of the External Libraries===========================
#=================================================================================
ExtLibDIR := Ext_Lib
# the name of the directory where to find all the external librairies (do not mind the ":")
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))
# first and foremost the presence of the directory containing the ext. libs. is tested as it could be done in a shell script

QD_DIR     = $(ExtLibDIR)/QDUtilLib
# the name and path where to find all the files of this specific external library (QDUtilLib). Rq: here its a lik-file pointing towards the QDUtilLib_loc file

OOPT       = 0
OOMP       = 1
LLAPACK    = 0
INT        = 4
# the parameters to compile the QDUtilLib with (cf makefile of QDUtilLib_loc)

ext_obj    = _$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
# the full name (with all the options) of the full static object file of the library (/!\ WITHOUT the .a extension and the prefix "libQD")
QDLIBA     = $(QD_DIR)/libQD$(ext_obj).a
# the full name and path (with all the extentions AND WITH the .a extention and the prefix "libQD" this time) : the actual static file of the QDUtilLib library
QDMOD_DIR  = $(QD_DIR)/OBJ/obj$(ext_obj)
# the name and path where to find all the .mod files of the QDUtilLib library

EXTLib     = $(QDLIBA)
EXTMod     = -I$(QDMOD_DIR)
# for the sake of legibility, the last two name/path (that is : the .a static lib. file and the directory where to find .mod files) with more intuitive names
# "-I" is for the gfortran compiler to understand that this is the path towards the .mod files for the compilation (gfortran argument syntax)


#=================================================================================
#==================Definition of other useful names in variables==================
#=========================3. For the gfortran compilation=========================
#=================================================================================
FFLAGS   = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp
# some useful optional arguments of the gfortran compilation command
FFLAGS   += -J$(MOD_DIR) $(EXTMod)
# "+=" works the same way as in python, C, etc
# to these optional aguments are added : "-J" which indicates the DIRECTORY (and not all the path of the .mod files) where to STORE the .mod files of the...
# ...lib. after compilation (at the contrary of "-I" which indicates where to FIND the needed ones); and the path towards the .mod files of the external...
# ...library(ies).
SRCFILES = Cavity_mode_m.f90 Operator_1D_m.f90 Operator_2D_m.f90 Total_hamiltonian_m.f90 Algebra_m.f90
# the list of all the .f90 source files OF THE LIBRARY (only the modules, not the test/app programs) to be compiled
OBJ0     = ${SRCFILES:.f90=.o}
# this syntax looks like a list slicing in python : change the .f90 string of the var to .o : it allows to change the extension of the file from .f90 to .o
# the list of all the names of the futur object files that will be created from the aforementioned library source files
OBJ      = $(addprefix $(OBJ_DIR)/, $(OBJ0))
# the command "addprefix <string to add as a prefix>, <string to which the prefix has to be added>" does exactly what it seems to
# the same list but, this time, with the complete path towards where the object files have to be stored

$(info ***********************************************************************)
$(info ***********************************************************************)
$(info ***********COMPILER:     $(FFC))
$(info ***********OBJ:          $(OBJ))
$(info ***********FFLAGS:       $(FFLAGS))
$(info ***********************************************************************)
$(info ***********************************************************************)
# the make is asked to display the chosen compiler, the list of the source/object files of the library to be created and the compilation arguments


#=================================================================================
#=========================Definition of the make commands=========================
#=======================1. the tests execution command "UT"=======================
#=================================================================================
# with make, a command (defined as hereinafter) is called from command-line with the syntax "make <name of the command>" or "make"
# here is defined the UT command that EXECUTE (and compile if necessary) all the programs testing the library as planned in the TESTS/ directory 
.PHONY: ut UT
# the ".PHONY <string1> <string2> <...>" make command indicates to make that the provided string are neither files nor directories and allows to use them...
# ... as key-words, ex: as command-line commands
UT ut: test_algebra.exe test_cavity_mode.exe #test_construct_op.exe test_action_op.exe
	./test_algebra.exe > $(OUTPUT_DIR)/test_algebra.log
	grep "Test" $(OUTPUT_DIR)/test_algebra.log
	./test_cavity_mode.exe < $(DATA_DIR)/data_tests.nml > $(OUTPUT_DIR)/test_cavity_mode.log
	grep "Test" $(OUTPUT_DIR)/test_cavity_mode.log
#	./test_construct_op.exe < $(DATA_DIR)/data_tests.nml > $(OUTPUT_DIR)/test_construct_op.log
#	grep "Test" $(OUTPUT_DIR)/test_construct_op.log
#	./test_action_op.exe < $(DATA_DIR)/data_tests.nml > $(OUTPUT_DIR)/test_action_op.log
#	grep "Test" $(OUTPUT_DIR)/test_action_op.log
	@echo "Done Tests"
# here is the actual definition of the command. It is declared as "UT" OR "ut" to make it cass-insensitive (both can be used to call it). NB: UT = Utility Test
# the first line instruction is understood by Make as "see these files". It will search the make file for where they are defined i.e. for their...
# ... dependancies, and create them as they are defined if they are too old. 
# the command will execute the executable files of the tests (and compile them first if they do not exist or if they are too old because that is the point...
# ... of using Make cf. definition of the files), using the default data_tests.nml data namelist (for now) and write the output in the .log file
# then it will display the output message of the test as written in the output file using the "grep" shell command. (So be careful programing your tests to...
# ... write the output massage as it is grepable)
# then if everything worked well it should print "Done Tests" using the "@echo" bash command.


#=================================================================================
#=========================Definition of the make commands=========================
#====================2. the application execution command "APP"===================
#=================================================================================
# this command will compile the library (create the .o and .mod files) and the tests (create the .o and .exe files) and create the static library .a file...
# ... BUT not execute anything !
.PHONY: app APP App
app APP App: $(MAIN).exe
	./$(MAIN).exe < $(DATA_DIR)/data_app.nml > $(OUTPUT_DIR)/$(MAIN).log


#=================================================================================
#=========================Definition of the make commands=========================
#===================3. the everything compilation command "all"===================
#=================================================================================
# this command will compile the library (create the .o and .mod files) and the tests (create the .o and .exe files) and create the static library .a file...
# ... BUT not execute anything !
.PHONY: all
all: $(LIBA) test_algebra.exe test_cavity_mode.exe test_construct_op.exe $(MAIN).exe
# Recall : LIBA = libMolecCav
# this instruction is understood by Make as "see these files". It will search the make file for where they are defined i.e. for their dependancies, and...
# ... create them as they are defined if they are too old. 


#=================================================================================
#=========================Definition of the make commands=========================
#===========4. the static library compilation and creation command "lib"==========
#=================================================================================
# this command will compile the library (create the .o and .mod files) into the static library file (create the .a file).
# the same command as "all" but without the compilation of the tests : the library only !  
.PHONY: lib
lib: $(LIBA)

$(LIBA): $(OBJ)
	ar -cr $(LIBA) $(OBJ)
	@echo "Done Library: "$(LIBA)
# the "ar -cr <name.a> <list of names of .o files>" command is the bash one to collate the provided object files into one static library file
# here is also added the definition of the static library file (instead of with the other file definitions below). cf. at the "Definition of the files"...
# ... section for more explanations


#=================================================================================
#=========================Definition of the make commands=========================
#====5. the static external libraries compilation and creation command "getlib"===
#=================================================================================
# this command will (install if possible and) compile the external libraries
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib
# with this instruction, make will actually execute the get_Lib.sh script as written in the directory of the corresponding external library, with the name...
# ... of the library (QDUtilLib) as an argument. This will install the library from GitHub 

$(QDLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "Done " $(QDLIBA)
# here is also added the definition of the static external library file (instead of with the other file definitions below). cf. at the "Definition of the...
# ... files" section for more explanations


#=================================================================================
#=========================Definition of the make commands=========================
#=================6. the cleaning commands "clean" and "cleanall"=================
#=================================================================================
# these commands are designed to get rid of all the files created at the compilation and execution of the library : object, modules, executables, statics,...
# ... outputs, etc
.PHONY: clean cleanall
clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f test*.exe
	rm -f $(MAIN).exe
	rm -f $(OUTPUT_DIR)/test_*.log
	rm -f $(OUTPUT_DIR)/$(MAIN).log
	@echo "Done cleaning objects, executables, and tests outputs"
# removes all the object files from the OBJ/ directory, all the executable files, and the output files from the tests

cleanall : clean
	rm -fr OBJ/obj*
	rm -f lib*.a
	cd Ext_Lib ; ./cleanlib
	@echo "Done all cleaning : objects, modules, statics, and same for external libraries"
# removes all the files from the OBJ/ directory (.o and .mod), all the static library files, and performs the cleaning the external libraries as defined...
#... in their makefile in addition to executing the "clean" command


#=================================================================================
#==================Definition of the files to be created by Make==================
#====================1. the main executable (the APP/ program)====================
#=================================================================================
$(MAIN).exe                    : $(OBJ_DIR)/$(MAIN).o $(LIBA)
	$(FFC) -o $(MAIN).exe $(FFLAGS) $(OBJ_DIR)/$(MAIN).o $(LIBA) $(EXTLib)
# this syntax define a file in Make : it shows its dependancies and then the compilation instructions
# the first line "<the file defined here> : <list of the files it depends on>" provides the dependancies. If Make has to treat this file, it will compare...
# ... its date of creation with the one of the files it depends on. If it has been created after them all nothing is done. Otherwise or if the file doesn't...
# ... exist, Make will create it using the instructions that follows
# the following lines are thus the creation instructions. They are usually compilation instruction but not necessarily
# /!\ they have to be indented with a TABULATION and not with spaces /!\.

$(OBJ_DIR)/$(MAIN).o           : $(MAIN_DIR)/$(MAIN).f90
	$(FFC) -c -o $(OBJ_DIR)/$(MAIN).o $(FFLAGS) $(MAIN_DIR)/$(MAIN).f90


#=================================================================================
#==================Definition of the files to be created by Make==================
#========================2. the tests (the TESTS/ programs)=======================
#=================================================================================
test_algebra.exe          : $(OBJ_DIR)/test_algebra.o $(LIBA)
	$(FFC) -o test_algebra.exe  $(FFLAGS) $(OBJ_DIR)/test_algebra.o $(LIBA) $(EXTLib)

$(OBJ_DIR)/test_algebra.o : $(TESTS_DIR)/test_algebra.f90
	$(FFC) -c -o $(OBJ_DIR)/test_algebra.o $(FFLAGS) $(TESTS_DIR)/test_algebra.f90

test_cavity_mode.exe           : $(OBJ_DIR)/test_cavity_mode.o $(LIBA)
	$(FFC) -o test_cavity_mode.exe  $(FFLAGS) $(OBJ_DIR)/test_cavity_mode.o $(LIBA) $(EXTLib)
# N.B. contrary to the next instruction, -o does mean here "link these object files into an executable one" (gfortran argument syntax)
# N.B. Question : je ne suis pas de pourquoi FFLAGS est nécéssaire pour le linking en .exe puisque je crois que les .mod ne sont pas utiles pour cette étape ?
# Recall: "FFLAGS also provides the path towards the .mod files : "-J$(MOD_DIR) $(EXTMod)"

$(OBJ_DIR)/test_cavity_mode.o  : $(TESTS_DIR)/test_cavity_mode.f90
	$(FFC) -c -o $(OBJ_DIR)/test_cavity_mode.o $(FFLAGS) $(TESTS_DIR)/test_cavity_mode.f90
# "-c" is the compilation and can take also as an argument "-o" which here does not mean "link into executable" but allows to choose the name of the thereby...
#... created object files. Here the path where they have to be stored (OBJ/obj) is added as a prefix to the name 

test_construct_op.exe          : $(OBJ_DIR)/test_construct_op.o $(LIBA)
	$(FFC) -o test_construct_op.exe  $(FFLAGS) $(OBJ_DIR)/test_construct_op.o $(LIBA) $(EXTLib)

$(OBJ_DIR)/test_construct_op.o : $(TESTS_DIR)/test_construct_op.f90
	$(FFC) -c -o $(OBJ_DIR)/test_construct_op.o $(FFLAGS) $(TESTS_DIR)/test_construct_op.f90

test_action_op.exe             : $(OBJ_DIR)/test_action_op.o $(LIBA)
	$(FFC) -o test_action_op.exe  $(FFLAGS) $(OBJ_DIR)/test_action_op.o $(LIBA) $(EXTLib)

$(OBJ_DIR)/test_action_op.o : $(TESTS_DIR)/test_action_op.f90
	$(FFC) -c -o $(OBJ_DIR)/test_action_op.o $(FFLAGS) $(TESTS_DIR)/test_action_op.f90

#=================================================================================
#==================Definition of the files to be created by Make==================
#=========================3. the library (the SRC modules)========================
#==============================3.1. The object files==============================
#=================================================================================
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FFC) $(FFLAGS) -o $@ -c $<

#$(OBJ_DIR)/Cavity_mode_m.o  : $(SRC_DIR)/Cavity_mode_m.f90
#	$(FFC) -c -o $(OBJ_DIR)/Cavity_mode_m.o $(FFLAGS) $(SRC_DIR)/Cavity_mode_m.f90

#$(OBJ_DIR)/Operator_1D_m.o  : $(SRC_DIR)/Operator_1D_m.f90
#	$(FFC) -c -o $(OBJ_DIR)/Operator_1D_m.o $(FFLAGS) $(SRC_DIR)/Operator_1D_m.f90

#$(OBJ_DIR)/Total_hamiltonian_m.o  : $(SRC_DIR)/Total_hamiltonian_m.f90
#	$(FFC) -c -o $(OBJ_DIR)/Total_hamiltonian_m.o $(FFLAGS) $(SRC_DIR)/Total_hamiltonian_m.f90

#$(OBJ_DIR)/Algebra_m.o  : $(SRC_DIR)/Algebra_m.f90
#	$(FFC) -c -o $(OBJ_DIR)/Algebra_m.o $(FFLAGS) $(SRC_DIR)/Algebra_m.f90

#=================================================================================
#==================Definition of the files to be created by Make==================
#=========================3. the library (the SRC modules)========================
#===========================3.2. The static library file==========================
#=================================================================================
# /!\ this file definition has been moved to the "3. the static library compilation and creation command "lib"" section /!\.


#=================================================================================
#==================Definition of the files to be created by Make==================
#=================4. the external libraries (the Ext_Lib modules)=================
#=================================================================================
# /!\ this file definition has been moved to the "4. the static external libraries compilation and creation command "getlib"" section /!\.


#=================================================================================
#=======================Definition of the files dependancies======================
#=================================================================================
# here are added the other dependancies between modules, that arre mandatory for a good compilation, but which do not lead to creation instructions. Just to...
# ... specify that one module needs another one
$(OBJ_DIR)/test_algebra.o        : $(LIBA)
$(OBJ_DIR)/test_cavity_mode.o    : $(LIBA)
$(OBJ_DIR)/test_construct_op.o   : $(LIBA)
$(OBJ_DIR)/test_action_op.o      : $(LIBA)

$(OBJ_DIR)/Operator_1D_m.o       : $(OBJ_DIR)/Algebra_m.o 
$(OBJ_DIR)/Operator_1D_m.o       : $(OBJ_DIR)/Cavity_mode_m.o 
$(OBJ_DIR)/Operator_2D_m.o       : $(OBJ_DIR)/Algebra_m.o 
$(OBJ_DIR)/Operator_2D_m.o       : $(OBJ_DIR)/Cavity_mode_m.o 
$(OBJ_DIR)/Operator_2D_m.o       : $(OBJ_DIR)/Operator_1D_m.o 
$(OBJ_DIR)/Total_hamiltonian_m.o : $(OBJ_DIR)/Cavity_mode_m.o 
$(OBJ_DIR)/Total_hamiltonian_m.o : $(OBJ_DIR)/Operator_1D_m.o 


$(OBJ_DIR)/$(MAIN).o             : $(LIBA)

$(OBJ)                           : | $(QDLIBA)
# the pipe means that only the EXISTENCE is tested and not the dates of the files. The logical test is True if $(QDLIBA) exists, EVEN if its creation date...
#... is more recent than the tested file

