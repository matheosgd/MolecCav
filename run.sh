#! /bin/bash

###############################################################################
###############################################################################
###############################################################################
## This script allows to build the library and run the tests in a similar m- ##
## anner than the makefile, despite not allowing as much flexibility. The F- ##
## ortran compiler cannot be changed and the date of creation of the files   ##
## are not compared to each other : all files are built each time. Contrary  ##
## to the make, only the "ut" command exists (not "UT") and so does the "ap- ##
## p" (not the "App" nor "APP").                                             ##
##                                                                           ##
## Syntax :                                                                  ##
## "./run.sh <arg1> <arg2> <arg3> <arg4>"; all the arguments being optionna- ##
## l, provided in any order.                                                 ##
## No more than four arguments are expected for now.                         ##
##                                                                           ##
## 1. One of the arguments is expected to be the name of the command to be   ##
## executed. The possible commands are :                                     ##
## 1.1. "ut" : if no other agruments are provided all the tests will be bui- ##
## lt and executed, and as a consequence, the library will be built too.     ##
## arg2 Can be provided. In that case, it is expected to be the name of a t- ##
## est WITHOUT THE EXTENSION, which will be the only test built and execute- ##
## d.                                                                        ## 
## 1.2. "app" : the application will be built and executed, and as a conseq- ##
## uence, the library will be built too. PLEASE do NOT change the name of t- ##
## he application source file !                                              ##
## 1.3. "all" : build the library (.mod; .o; .a files) and the executable f- ##
## ile of all the tests but do not execute anything.                         ##
## 1.4. "lib" : same command as "all" but without the test .i.e. build the   ##
## library.                                                                  ##
## 1.5. "getlib" :  with this instruction, the script will actually execute  ##
## the get_Lib.sh script as written in the directory of the corresponding e- ##
## xternal library, with the name of the library (QDUtilLib) as an argument. ##
## This will install the library from GitHub.                                ##
## 1.6. "clean" : removes all the object files from the OBJ/ directory, all  ##
## the executable files, and the output files from the tests.                ##
## 1.7. "cleanall" : removes all the files from the OBJ/ directory (.o and   ##
# .mod), all the static library files, and performs the cleaning of the ext- ##
## ernal libraries as defined in their makefile.                             ##
## 1.8. The default is arg1 = "ut"; arg2 = "".                               ##
##                                                                           ##
## 2. The other arguments are expected to be :                               ##
## 2.1. either the name /!\ without the extension /!\ of the test to be exe- ##
## cuted (in the 1.1. case) /!\ in the "test_<test_name>" format /!\...      ##
## 2.2. ...or the name /!\ without the extension /!\ of the data file (the   ##
## namelist) to be read (case 1.1., 1.2., 1.3.) /!\ in the "data_<data_name>"##
## format /!\...                                                             ##
## 2.3. or the name of the obj/ subdirectory to be cleaned or filled /!\ in  ##
## the "OBJ/<name_of_the_subdirectory>" format.                              ##
###############################################################################
###############################################################################
###############################################################################

echo "################################################################################"
echo "################################################################################"


#------------------------------------------------------------------------------
#---------------------------------0. Utilities---------------------------------
#---------------------------0.1 Useful variable names--------------------------
#------------------------------------------------------------------------------
FFC="gfortran"                                                                 # the fortran compiler
FFLAGS="-Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp"
                                                                               # some useful options for the compiler

MODULES_LIB=("Cavity_mode_m" "Operator_1D_m" "Total_hamiltonian_m" "Algebra_m" "Operator_2D_m")
SRC_LIB=(${MODULES_LIB[@]/%/.f90})
OBJ_LIB=(${MODULES_LIB[@]/%/.o})

LIB="libMolecCav"
LIBA="$LIB.a"

TESTS=("test_algebra" "test_cavity_mode" "test_construct_op_1D" "test_action_op_1D" "test_action_total_H_1p1D" "test_construct_total_H_1p1D")
SRC_TESTS=(${TESTS[@]/%/.f90})                                                 # parenthesis must be added here to make bash know that the next var is also an array
OBJ_TESTS=(${TESTS[@]/%/.o})                                                   
EXE_TESTS=(${TESTS[@]/%/.exe})  

MAIN_DIR="APP"
MAIN="App_MolecCav"
SRC_MAIN="$MAIN.f90"
OBJ_MAIN="$MAIN.o"
EXE_MAIN="$MAIN.exe"
OUT_MAIN="$MAIN.log"

OUTPUT_DIR="OUT"

DATA_DIR="DATA"

ExtLibDIR="Ext_Lib"
OOPT="0"
OOMP="1"
LLAPACK="0"
INT="4"
EXTMod="-I$ExtLibDIR/QDUtilLib/OBJ/obj_${FFC}_opt${OOPT}_omp${OOMP}_lapack${LLAPACK}_int${INT}"
QDLIBA="$ExtLibDIR/QDUtilLib/libQD_${FFC}_opt${OOPT}_omp${OOMP}_lapack${LLAPACK}_int${INT}.a"
EXTLib="$QDLIBA"


#------------------------------------------------------------------------------
#-----------------------0.2 Declaration of the functions-----------------------
#------------------------------------------------------------------------------
Claim()
{
  if [ ! "$1" = "/end" ]
  then
    Nb_characters="${#1}"
    blank_line="                                                                                "

    Nb_lines=$((1+$Nb_characters/73))

    for ((i=0 ; $((Nb_lines-1)) - $i ; i++))
    do
      line="## ${1:$((i*73)):73}- ##"
      echo -e "$line"
    done

    remaining_characters="$(($Nb_characters-73*$Nb_lines+73))"
    Nb_spaces="$((74-$remaining_characters))"
    Demi_Nb_spaces="$(($Nb_spaces/2))"
    remains="$(($Nb_spaces%2))"
  
    line="## ${blank_line:0:$((Demi_Nb_spaces+$remains))}${1:$((73*Nb_lines-73))}${blank_line:0:$Demi_Nb_spaces} ##"
    echo -e "$line"
  fi 

  if [ "$1" = "/end" -o "$2" = "/end" ]
  then
      echo "##----------------------------------------------------------------------------##"
  fi
}

Build_QDLIBA()
{
  cd ~/MolecCav
  Claim "Running the ./get_Lib.sh script of QDUtilLib :"
  Claim ":"
  Claim ":"
	cd $ExtLibDIR ; ./get_Lib.sh QDUtilLib
	cd QDUtilLib
  make lib FC=$FFC OPT=$OOPT OMP=$OOMP LAPACK=$LLAPACK INT=$INT ExtLibDIR=$ExtLibDIR CompilersDIR=$CompilersDIR
  cd ~/MolecCav/
	test -f $QDLIBA || (echo $QDLIBA "does not exist" ; exit 1)
  Claim ":"
  Claim ":"
  Claim "Done $QDLIBA" "/end"
}

Build_obj_lib()
{
  cd ~/MolecCav
 	$FFC -c -o ${OBJ[0]} $FFLAGS ${SRC_FILES[0]}                                           # OBJ[0] = OBJ/obj/Cavity_mode_m.o 
                                                                                         # SRC_FILES[0] = OBJ/obj/Cavity_mode_m.f90
  $FFC -c -o ${OBJ[3]} $FFLAGS ${SRC_FILES[3]}                                           # OBJ[2] = OBJ/obj/Algebra_m.o
                                                                                         # SRC_FILES[2] = OBJ/obj/Algebra_m.f90    |/!\ Algebra has to be compiled         /!\
  $FFC -c -o ${OBJ[1]} $FFLAGS ${SRC_FILES[1]}                                           # OBJ[1] = OBJ/obj/Cavity_mode_m.o        |/!\ before operator_1d (used in op_1d) /!\
                                                                                         # SRC_FILES[1] = OBJ/obj/Cavity_mode_m.f90
  $FFC -c -o ${OBJ[4]} $FFLAGS ${SRC_FILES[4]}                                           # OBJ[4] = OBJ/obj/Algebra.o
                                                                                         # SRC_FILES[2] = OBJ/obj/Algebra.f90
  $FFC -c -o ${OBJ[2]} $FFLAGS ${SRC_FILES[2]}                                           # OBJ[2] = OBJ/obj/Total_hamiltonian_m.o
                                                                                         # SRC_FILES[2] = OBJ/obj/Total_hamiltonian_m.f90

  for file in ${OBJ[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"
}

Build_liba()
{
  ar -cr $LIBA ${OBJ[@]}
  Claim "Done Library $LIBA ..."
  Claim "...using ${OBJ[*]}" "/end"
}

Build_lib()
{
  Build_QDLIBA
  Build_obj_lib
  Build_liba
}

Build_tests()
{
	$FFC -c -o ${TESTS_OBJ_FILES[0]} $FFLAGS ${TESTS_SRC_FILES[0]}               #it turned out -JOBJ/obj is enough and do not need -IOBJ/obj in addition
	$FFC -o ${EXE_TESTS[0]}  $FFLAGS ${TESTS_OBJ_FILES[0]} $LIBA $EXTLib
	$FFC -c -o ${TESTS_OBJ_FILES[1]} $FFLAGS ${TESTS_SRC_FILES[1]}
	$FFC -o ${EXE_TESTS[1]}  $FFLAGS ${TESTS_OBJ_FILES[1]} $LIBA $EXTLib
	$FFC -c -o ${TESTS_OBJ_FILES[2]} $FFLAGS ${TESTS_SRC_FILES[2]}
	$FFC -o ${EXE_TESTS[2]}  $FFLAGS ${TESTS_OBJ_FILES[2]} $LIBA $EXTLib
	$FFC -c -o ${TESTS_OBJ_FILES[3]} $FFLAGS ${TESTS_SRC_FILES[3]}
	$FFC -o ${EXE_TESTS[3]}  $FFLAGS ${TESTS_OBJ_FILES[3]} $LIBA $EXTLib
	$FFC -c -o ${TESTS_OBJ_FILES[4]} $FFLAGS ${TESTS_SRC_FILES[4]}
	$FFC -o ${EXE_TESTS[4]}  $FFLAGS ${TESTS_OBJ_FILES[4]} $LIBA $EXTLib
	$FFC -c -o ${TESTS_OBJ_FILES[5]} $FFLAGS ${TESTS_SRC_FILES[5]}
	$FFC -o ${EXE_TESTS[5]}  $FFLAGS ${TESTS_OBJ_FILES[5]} $LIBA $EXTLib

  for file in ${TESTS_OBJ_FILES[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"

  for file in ${EXE_TESTS[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"
}

Build_app()
{
  $FFC -c -o $OBJ_DIR/$OBJ_MAIN $FFLAGS $MAIN_DIR/$SRC_MAIN
  Claim "Done $OBJ_MAIN" "/end"
  $FFC -o $EXE_MAIN $FFLAGS $OBJ_DIR/$OBJ_MAIN $LIBA $EXTLib                   # from the most internal module to the least internal one
  Claim "Done $EXE_MAIN" "/end"
}


#------------------------------------------------------------------------------
#-------------------------1. Recovery of the arguments-------------------------
#--------------------------1.0. The default parameters-------------------------
#------------------------------------------------------------------------------
command="ut"
data_file="data_tests"
test_name=""
OBJ_DIR="OBJ/obj"
MOD_DIR="$OBJ_DIR"
SRC_DIR="SRC"
TESTS_DIR="TESTS"

FFLAGS+=" -J$MOD_DIR $EXTMod"

if [ -z "$1" ]                                             # z=!n; tests if the value is an empty string.
then
  Claim "No instruction provided !"
  Claim "The script will be run with the default parameters"
  Claim "cf. \"Syntaxe\""
  Claim "/end"
fi


#------------------------------------------------------------------------------
#-------------------------1. Recovery of the arguments-------------------------
#-------------------------------1.1. The command-------------------------------
#------------------------------------------------------------------------------
case "$1" in
"ut" | "all" | "lib" | "getlib" | "clean" | "cleanall" | "app")
  command="$1";;
*)
  Claim "The first argument is not a command"
esac
case "$2" in
"ut" | "all" | "lib" | "getlib" | "clean" | "cleanall" | "app")
  command="$2";;
*)
  Claim "The second argument is not a command"
esac
case "$3" in
"ut" | "all" | "lib" | "getlib" | "clean" | "cleanall" | "app")
  command="$3";;
*)
  Claim "The third argument is not a command"
esac
case "$4" in
"ut" | "all" | "lib" | "getlib" | "clean" | "cleanall" | "app")
  command="$4";;
*)
  Claim "The fourth argument is not a command";;
esac
  Claim "/end"
#echo "\$1=$1" #echo "\$command=$command"


#------------------------------------------------------------------------------
#------------------------------1.2. The test name-----------------------------
#------------------------------------------------------------------------------
case "test_" in
"${1:0:5}")
  test_name="$1.f90";;
"${2:0:5}")
  test_name="$2.f90";;
"${3:0:5}")
  test_name="$3.f90";;
"${4:0:5}")
  test_name="$4.f90";;
*) 
  if [ $command = "ut" ]
  then
    Claim "No test name provided !"
    Claim "All of them will be executed" "/end"
  else
    Claim "No test name provided !" "/end"
  fi;;
esac


#------------------------------------------------------------------------------
#-------------------------------1.3. The namelist------------------------------
#------------------------------------------------------------------------------
case "data_" in
"${1:0:5}")
  data_file="$1.nml";;
"${2:0:5}")
  data_file="$2.nml";;
"${3:0:5}")
  data_file="$3.nml";;
"${4:0:5}")
  data_file="$4.nml";;
*) 
  case "$command" in
  "lib" | "all" | "ut" | "app")
    Claim "No data_file provided !"
    Claim "The default data_tests.nml will be used" "/end";;
  "clean" | "cleanall" | "getlib")
    Claim "No data_file provided !" "/end";;
  *) 
    Claim "something is weird 1"
    echo "################################################################################"
    echo "################################################################################"
    exit 1;;
  esac;;
esac


#------------------------------------------------------------------------------
#-------------------------------1.4. The OBJ_DIR-------------------------------
#------------------------------------------------------------------------------
case "OBJ/" in
"${1:0:4}")
  OBJ_DIR="$1"
  MOD_DIR="$OBJ_DIR";;
"${2:0:4}")
  OBJ_DIR="$2"
  MOD_DIR="$OBJ_DIR";;
"${3:0:4}")
  OBJ_DIR="$3"
  MOD_DIR="$OBJ_DIR";;
"${4:0:4}")
  OBJ_DIR="$4"
  MOD_DIR="$OBJ_DIR";;
*)
  if [ ! $command = "getlib" ]
  then
    Claim "No OBJ/subdirectory provided !"
    Claim "The \"obj\" default one will be used" "/end"
  else
    Claim "No OBJ/subdirectory provided !" "/end"
  fi;;
esac

if [ ! -d "$OBJ_DIR" ]                                                             # d teste l'existence du directory "<...>"
then
  mkdir "$OBJ_DIR"
  Claim "$OBJ_DIR directory was not here yet : created" "/end"
else
  Claim "$OBJ_DIR directory already created : ok" "/end"
fi


#------------------------------------------------------------------------------
#----------------------------------1.5. Sum-Up---------------------------------
#------------------------------------------------------------------------------
Claim "SUM UP :"
case "$command" in
"lib")
  Claim "This script will run the command .................... $command ..."
  Claim "...using the data of the namelist .......... $data_file.nml ..."
  Claim "...and using, for the .o and .mod files, the subdirectory ..... $OBJ_DIR...";;
"all")
  Claim "This script will run the command .................... $command ..."
  Claim "...using the data of the namelist .......... $data_file.nml ..."
  Claim "...and using, for the .o and .mod files, the subdirectory ..... $OBJ_DIR...";;
"ut")
  Claim "This script will run the command .................... $command ..."
  Claim "...using the data of the namelist .......... $data_file.nml ..."
  Claim "...using, for the .o and .mod files, the directory ..... $OBJ_DIR ..."
  if [ -n "$test_name" ]
  then
    Claim "...and executing only the test ............... $test_name"
  else
  Claim "...and executing all the tests."
  fi;;
"app")
  Claim "This script will run the command .................... $command ..."
  Claim "...using the data of the namelist .......... data_app.nml ..."
  Claim "...using, for the .o and .mod files, the directory ..... $OBJ_DIR ...";;
"getlib")
  Claim "This script will run the command .................... $command...";;
"clean")
  Claim "This script will run the command .................... $command...";;
"cleanall")
  Claim "This script will run the command .................... $command...";;
*) 
  Claim "something is weird 2"
  Claim "################################################################################"
  Claim "################################################################################"
  exit 1;;
esac
  Claim "/end"


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#------------------------------------------------------------------------------
OBJ=(${OBJ_LIB[@]/#/$OBJ_DIR/})
SRC_FILES=(${SRC_LIB[@]/#/$SRC_DIR/})
TESTS_SRC_FILES=(${SRC_TESTS[@]/#/$TESTS_DIR/})
TESTS_OBJ_FILES=(${OBJ_TESTS[@]/#/$OBJ_DIR/})
MAIN_OBJ_FILES=$OBJ_DIR/$OBJ_MAIN
MAIN_SRC_FILES=$MAIN_DIR/$MAIN
#echo "\${OBJ[@]} = ${OBJ[@]}"
#echo "\${SRC_FILES[@] = ${SRC_FILES[@]}"
#echo "TESTS_FILES = ${TESTS_SRC_FILES[@]}"
#echo "OBJ_FILES = ${TESTS_OBJ_FILES[@]}"

case "$command" in


#------------------------------------------------------------------------------
#----------------------------2.1. The clean command----------------------------
#------------------------------------------------------------------------------
"clean") 
	rm -f $OBJ_DIR/*.o
	rm -f test*.exe
  rm -f $EXE_MAIN
	rm -f $OUTPUT_DIR/test_*.log
	rm -f $OUTPUT_DIR/$OUT_MAIN
	Claim "Done cleaning objects, executables, and tests outputs"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#---------------------------2.1. The cleanall command--------------------------
#------------------------------------------------------------------------------
"cleanall")
  ./run.sh clean 	
  rm -fr OBJ/obj*
	rm -f lib*.a
	cd Ext_Lib ; ./cleanlib
	Claim "Done all cleaning :"
  Claim "objects, modules, statics, and same for external libraries"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#----------------------------2.1. The getlib command---------------------------
#------------------------------------------------------------------------------
"getlib")
 	cd $ExtLibDIR
  ./get_Lib.sh QDUtilLib
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#-----------------------------2.1. The lib command-----------------------------
#------------------------------------------------------------------------------
"lib")
  Build_lib
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#-----------------------------2.1. The all command-----------------------------
#------------------------------------------------------------------------------
"all")
  Build_lib
  #echo "${TESTS_OBJ_FILES[@]} ${TESTS_SRC_FILES[@]} ${EXE_TESTS[@]}"
  Build_tests
  Build_app
  Claim "Done library, tests executables, and application executable"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#------------------------------2.1. The ut command-----------------------------
#------------------------------------------------------------------------------
"ut")
  if [ ! -d "$OUTPUT_DIR" ]                                                    # d teste l'existence du directory "<...>"
  then
    mkdir "$OUTPUT_DIR"
    Claim "$OUTPUT_DIR/ directory was not here yet : created" "/end"
  else
    Claim "$OUTPUT_DIR directory already created : ok" "/end"
  fi

  Build_lib
  Build_tests
  #./$test_name.exe < $data_file.nml > $data_file.out
  ./test_algebra.exe > "$OUTPUT_DIR/"test_algebra.log
  Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_file_sca_pdt.log)"
  Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_file_norm.log)"
  Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_file_normalztn.log)"
	./test_cavity_mode.exe < ${DATA_DIR}/data_tests.nml > "$OUTPUT_DIR/"test_cavity_mode.log
	Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_cavity_mode.log)"
	./test_construct_op_1D.exe < ${DATA_DIR}/data_tests.nml > "$OUTPUT_DIR/"test_construct_op_1D.log
	Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_construct_op_1D.log)"
	./test_action_op_1D.exe < ${DATA_DIR}/data_tests.nml > "$OUTPUT_DIR/"test_action_op_1D.log
	Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_action_op_1D.log)"
	./test_action_total_H_1p1D.exe < ${DATA_DIR}/data_tests.nml > "$OUTPUT_DIR/"test_action_total_H_1p1D.log
	Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_action_total_H_1p1D.log)"
	./test_construct_total_H_1p1D.exe < ${DATA_DIR}/data_tests.nml > "$OUTPUT_DIR/"test_construct_total_H_1p1D.log
	Claim "$(grep "Number of error(s)" "$OUTPUT_DIR/"test_construct_total_H_1p1D.log)"
	Claim "Done Tests"

  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#------------------------------2.1. The ut command-----------------------------
#------------------------------------------------------------------------------
"app")
  if [ ! -d "$OUTPUT_DIR" ]                                                    # d teste l'existence du directory "<...>"
  then
    mkdir "$OUTPUT_DIR"
    Claim "$OUTPUT_DIR/ directory was not here yet : created" "/end"
  else
    Claim "$OUTPUT_DIR directory already created : ok" "/end"
  fi

  Build_lib
  Build_app
	./$EXE_MAIN < ${DATA_DIR}/data_app.nml > "$OUTPUT_DIR/"$OUT_MAIN
	Claim "Done application"

  echo "################################################################################"
  echo "################################################################################"
  exit 0;;

*)
  Claim "something is weird 3"
  echo "################################################################################"
  echo "################################################################################"
  exit 1;;
esac

