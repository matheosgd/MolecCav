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
# UPPERCASE VAR = INTERNAL TO THE SCRIPT 
# LOWERCASE VAR = RECOVERED FROM ARGUMENTS (except for OBJ_DIR if changed)

#------------------------------------------------------------------------------
#---------------------------------0. Utilities---------------------------------
#---------------------------0.1 Useful variable names--------------------------
#------------------------------------------------------------------------------
FFC="gfortran"                                                                 # the fortran compiler

MODULES_DIR="SRC"
MODULES=("P1D_basis_m" "P1D_lanczos_m" "P1D_operators_m" "P1D_perturbations_m" "P1D_propagation_m" "P1D_wavepacket_m")
MODULES_SRC=(${MODULES[@]/%/.f90})                                             # /%/ add the following suffix to the all (because of "@") the elements of the targeted tabular

LIB="libPropagation_1D"
LIBA="$LIB.a"

TESTS_DIR="TESTS"
TESTS=("nope_there_is_none_up_to_now")
TESTS_SRC=(${TESTS[@]/%/.f90})                                                 # parenthesis must be added here to make bash know that the nex var is also an array
TESTS_EXE=(${TESTS[@]/%/.exe})  

MAIN_DIR="APP"
MAIN=("Main" "App_Propagation_1D")
MAIN_SRC=(${MAIN[@]/%/.f90})                                                   # parenthesis must be added here to make bash know that the nex var is also an array
MAIN_EXE=(${MAIN[@]/%/.exe})  

OBJ_DIR="OBJ/obj"
MOD_DIR="$OBJ_DIR"
MODULES_OBJ=(${MODULES[@]/%/.o})
TESTS_OBJ=(${TESTS[@]/%/.o})                                                   
MAIN_OBJ=(${MAIN[@]/%/.o})                                                   

DATA_DIR="DATA"

OUTPUT_DIR="OUT"
TESTS_OUT=(${TESTS[@]/%/.log})  
MAIN_OUT=(${MAIN[@]/%/.log})  

ExtLibDIR="Ext_Lib"
OOPT="0"
OOMP="1"
LLAPACK="0"
INT="4"
EXTMod="-I$ExtLibDIR/QDUtilLib/OBJ/obj_${FFC}_opt${OOPT}_omp${OOMP}_lapack${LLAPACK}_int${INT}"
QDLIBA="$ExtLibDIR/QDUtilLib/libQD_${FFC}_opt${OOPT}_omp${OOMP}_lapack${LLAPACK}_int${INT}.a"
EXTLib="$QDLIBA"

FFLAGS="-Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp" # some useful options for the compiler
FFLAGS+=" -J$MOD_DIR $EXTMod"


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
      echo "$line"
    done

    remaining_characters="$(($Nb_characters-73*$Nb_lines+73))"
    Nb_spaces="$((74-$remaining_characters))"
    Demi_Nb_spaces="$(($Nb_spaces/2))"
    remains="$(($Nb_spaces%2))"
  
    line="## ${blank_line:0:$((Demi_Nb_spaces+$remains))}${1:$((73*Nb_lines-73))}${blank_line:0:$Demi_Nb_spaces} ##"
    echo "$line"
  fi 

  if [ "$1" = "/end" -o "$2" = "/end" ]
  then
      echo "##----------------------------------------------------------------------------##"
  fi
}

Build_QDLIBA()
{
  cd ~/stageicp/Propagation_1D/
  Claim "Running the ./get_Lib.sh script of QDUtilLib :"
  Claim ":"
  Claim ":"
	cd "$ExtLibDIR" ; ./get_Lib.sh QDUtilLib
	cd QDUtilLib
  Claim ":"
  Claim ":"
  Claim "Running the makefile of QDUtilLib to compile the library :"
  Claim ":"
  Claim ":"
    make lib FC=$FFC OPT=$OOPT OMP=$OOMP LAPACK=$LLAPACK INT=$INT ExtLibDIR=$ExtLibDIR CompilersDIR=$CompilersDIR
  Claim ":"
  Claim ":"
  cd ~/stageicp/Propagation_1D/
  test -f $QDLIBA || (echo $QDLIBA "does not exist" ; exit 1)
  Claim "Done $QDLIBA" "/end"
}

Build_MODULES_OBJ()
{
  cd ~/stageicp/Propagation_1D/                                                # /!\ THE ORDER IN WHICH THE MODULES ARE COMPILED IS CRUCIAL FOR POTENTIAL INTERCONNECTIONS /!\
  $FFC -c -o $OBJ_DIR/${MODULES_OBJ[0]} $FFLAGS $MODULES_DIR/${MODULES_SRC[0]} # MODULES_OBJ[0] = P1D_basis_m.o
                                                                               # MODULES_SRC[0] = P1D_basis_m.f90
  $FFC -c -o $OBJ_DIR/${MODULES_OBJ[2]} $FFLAGS $MODULES_DIR/${MODULES_SRC[2]}
  $FFC -c -o $OBJ_DIR/${MODULES_OBJ[5]} $FFLAGS $MODULES_DIR/${MODULES_SRC[5]}
  $FFC -c -o $OBJ_DIR/${MODULES_OBJ[1]} $FFLAGS $MODULES_DIR/${MODULES_SRC[1]}
  $FFC -c -o $OBJ_DIR/${MODULES_OBJ[4]} $FFLAGS $MODULES_DIR/${MODULES_SRC[4]}
  $FFC -c -o $OBJ_DIR/${MODULES_OBJ[3]} $FFLAGS $MODULES_DIR/${MODULES_SRC[3]}

  for file in ${MODULES_OBJ[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"
}

Build_liba()
{
  ar -cr $LIBA ${MODULES_OBJ[@]/#/${OBJ_DIR}/}                                 # /#/ add the following prefix to the all (because of "@") the elements of the targeted tabular
  Claim "Done Library $LIBA ..."
  Claim "...using ${MODULES_OBJ[*]}" "/end"
}

Build_lib()
{
  Build_QDLIBA
  Build_MODULES_OBJ
  Build_liba
}

Build_tests()
{
  $FFC -c -o $OBJ_DIR/${TESTS_OBJ[0]} $FFLAGS $TESTS_DIR/${TESTS_SRC[0]}       #it turned out -JOBJ/obj is enough and do not need -IOBJ/obj in addition

  for file in ${TESTS_OBJ[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"

  $FFC -o ${TESTS_EXE[0]} $FFLAGS $OBJ_DIR/${TESTS_OBJ[0]} $LIBA $EXTLib

  for file in ${TESTS_EXE[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"
}

Build_app()
{
  $FFC -c -o $OBJ_DIR/${MAIN_OBJ[0]} $FFLAGS $MAIN_DIR/${MAIN_SRC[0]}       #it turned out -JOBJ/obj is enough and do not need -IOBJ/obj in addition
  $FFC -c -o $OBJ_DIR/${MAIN_OBJ[1]} $FFLAGS $MAIN_DIR/${MAIN_SRC[1]}

  for file in ${MAIN_OBJ[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"

  $FFC -o ${MAIN_EXE[0]} $FFLAGS $OBJ_DIR/${MAIN_OBJ[0]} $LIBA $EXTLib
  $FFC -o ${MAIN_EXE[1]} $FFLAGS $OBJ_DIR/${MAIN_OBJ[1]} $LIBA $EXTLib

  for file in ${MAIN_EXE[@]}
  do
    Claim "Done $file"
  done
  Claim "/end"
}


#------------------------------------------------------------------------------
#-------------------------1. Recovery of the arguments-------------------------
#--------------------------1.0. The default parameters-------------------------
#------------------------------------------------------------------------------
command="ut"
test_name=""                                                                   # /!\ useful only for the execution i.e. whithin the "ut" command
app_name=""                                                                    # /!\ useful only for the execution i.e. whithin the "app" command

if [ -z "$1" ]                                                                 # z=!n; tests if the value is an empty string.
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
"getlib" | "lib" | "all" | "ut" | "app" | "clean" | "cleanall")
  command="$1";;
*)
  Claim "The first argument is not a command"
esac
case "$2" in
"getlib" | "lib" | "all" | "ut" | "app" | "clean" | "cleanall")
  command="$2";;
*)
  Claim "The second argument is not a command"
esac
case "$3" in
"getlib" | "lib" | "all" | "ut" | "app" | "clean" | "cleanall")
  command="$3";;
*)
  Claim "The third argument is not a command"
esac
case "$4" in
"getlib" | "lib" | "all" | "ut" | "app" | "clean" | "cleanall")
  command="$4";;
*)
  Claim "The fourth argument is not a command";;
esac
  Claim "/end"
#echo "\$1=$1" #echo "\$command=$command"


#------------------------------------------------------------------------------
#-------------------------------1.2. The OBJ_DIR-------------------------------
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
#----------------------------1.3. The data namelist----------------------------
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
  "ut")
    Claim "No data_file provided !"
    Claim "The default data_tests.nml will be used" "/end"
    data_file="data_tests.nml";;
  "app")
    Claim "No data_file provided !"
    Claim "The default data_app.nml will be used"
    Claim "WARNING if your intent is to execute your Main !" "/end"
    data_file="data_app.nml";;
  "getlib" | "lib" | "all" | "clean" | "cleanall")
    Claim "No data_file provided, but no need anyway for this command" "/end";;
  *) 
    Claim "something is weird : either the \"data_file\" name or the \"command\" is not recognized"
    echo "################################################################################"
    echo "################################################################################"
    exit 1;;
  esac;;
esac


#------------------------------------------------------------------------------
#------------------------------1.4. The test name-----------------------------
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
    Claim "No test name provided, but this command won't execute any test anyway" "/end"
  fi;;
esac


#------------------------------------------------------------------------------
#-------------------------------1.5. The app name------------------------------
#------------------------------------------------------------------------------
case "Main" in
"$1" | "$2" | "$3" | "$4")
  app_name="Main";;
*) 
  if [ $command = "app" ]
  then
    Claim "The name \"Main\" has not been provided for the app file to be executed"
    Claim "If no other app name has been provided, the default App_Propagation_1D will be"
    Claim ""
  else
    Claim "The name \"Main\" has not been provided, but this command won't execute any app anyway"
    Claim ""
  fi;;
esac
case "App_Propagation_1D" in
"$1" | "$2" | "$3" | "$4")
  app_name="App_Propagation_1D";;
*) 
  if [ $command = "app" ]
  then
    Claim "The name \"App_Propagation_1D\" has not been provided for the app file to be executed"
    Claim "If no other app name has been provided, the default App_Propagation_1D will be"
  else
    Claim "The name \"App_Propagation_1D\" has not been provided, but this command won't execute any app anyway"
  fi;;
esac
if [ -z "$app_name" -a $command = "app" ]
then
    Claim "No app name has been provided for the app file to be executed"
    Claim "The default App_Propagation_1D will be" "/end"
    app_name="App_Propagation_1D"
else
    Claim "No app name has been provided, but this command won't execute any app anyway" "/end"
fi


#------------------------------------------------------------------------------
#----------------------------------1.5. Sum-Up---------------------------------
#------------------------------------------------------------------------------
Claim "SUM UP :"
case "$command" in
"lib")
  Claim "This script will run the command .................... $command ..."
  #Claim "...compiling the modules ${MODULES[*]} ..."
  Claim "...and using, for the .o and .mod files, the subdirectory ..... $OBJ_DIR";;
"all")
  Claim "This script will run the command .................... $command ..."
  Claim "...and using, for the .o and .mod files, the subdirectory ..... $OBJ_DIR";;
"ut")
  Claim "This script will run the command .................... $command ..."
  Claim "...using the data of the namelist .......... $data_file ..."
  Claim "...using, for the .o and .mod files, the directory ..... $OBJ_DIR ..."
  if [ -n "$test_name" ]                                                       #-n is true if the following variable is not an empty string
  then
    Claim "...and executing only the test ............... $test_name"
  else
  Claim "...and executing all the tests."
  fi;;
"app")
  Claim "This script will run the command .................... $command ..."
  Claim "...using the data of the namelist ............. $data_file ..."
  Claim "...using, for the .o and .mod files, the directory ..... $OBJ_DIR";;
"getlib")
  Claim "This script will run the command .................... $command";;
"clean")
  Claim "This script will run the command .................... $command..."
  Claim "...targeting the subdirectory .................... $OBJ_DIR";;
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
#OBJ=(${OBJ_LIB[@]/#/$OBJ_DIR/})
#SRC_FILES=(${SRC_LIB[@]/#/$SRC_DIR/})
#TESTS_SRC_FILES=(${SRC_TESTS[@]/#/$TESTS_DIR/})
#TESTS_OBJ_FILES=(${OBJ_TESTS[@]/#/$OBJ_DIR/})
#MAIN_OBJ_FILES=$OBJ_DIR/$OBJ_MAIN
#MAIN_SRC_FILES=$MAIN_DIR/$MAIN
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
  rm -f ${TESTS_EXE[*]}
  rm -f ${MAIN_EXE[*]}
  rm -f ${TESTS_OUT[@]/#/${OUTPUT_DIR}/}
  #echo "$OUTPUT_DIR;" "$OUTPUT_DIR/${MAIN_OUT[*]}" ${TESTS_OUT[@]/#/${OUTPUT_DIR}/}
  INTERMEDIARY=(${MAIN[@]/#/${OUTPUT_DIR}/*})                                  # we want here to add as a prefix "OUT/*"... from $MAIN and not MAIN_OUT because .log would bother to add the "*" suffix
  #echo "${INTERMEDIARY[*]}" #echo "${INTERMEDIARY[@]/%/*.log}"
  rm -f ${INTERMEDIARY[@]/%/*.log}                                             # here the suffix "*.log" is eventually added. A intermediary variable is mandatory since adding a prefix AND a suffix must be done in two steps in shell/bash. N.B. the ".log" is not mandatory but safe
  Claim "Done cleaning objects, executables, and outputs"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#---------------------------2.1. The cleanall command--------------------------
#------------------------------------------------------------------------------
"cleanall")
  ./run.sh clean 	
  rm -fr OBJ/*
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
  Claim "Running the get_Lib.sh script from QDUtilLib :"
  Claim ":"
  Claim ":"
    ./get_Lib.sh QDUtilLib
  Claim ":"
  Claim ":"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#-----------------------------2.1. The lib command-----------------------------
#------------------------------------------------------------------------------
"lib")
  Build_lib
  Claim "Done library $LIBA"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#-----------------------------2.1. The all command-----------------------------
#------------------------------------------------------------------------------
"all")
  Build_lib
  #Build_tests
  Build_app
  Claim "Done library $LIBA..."
  Claim "...tests executables ${TESTS_EXE[*]}..."
  Claim "...and applications executables ${MAIN_EXE}"
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
#######################################FOR NOW THERE IS NO TESTS TO COMPILE###################################
exit 0 # when it will be tests just remove this exit !

  Build_tests 

  if [ -n "$test_name" ]
  then
    ./${test_name}.exe < $DATA_DIR/$data_file > "$OUTPUT_DIR/"${test_name}.log           #N.B. the "{}" are not mandatory here
  else
    ./${TESTS_EXE[0]} < $DATA_DIR/$data_file > "$OUTPUT_DIR/"${TESTS_OUT[0]}
  fi
  Claim "$(grep "Test" "$OUTPUT_DIR/"${TESTS_OUT[0]})"
  
  Claim "Done Tests"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;
#######################################FOR NOW THERE IS NO TESTS TO COMPILE###################################


#------------------------------------------------------------------------------
#--------------------------2. Execution of the command-------------------------
#------------------------------2.1. The app command----------------------------
#------------------------------------------------------------------------------
"app")
  if [ ! -d "$OUTPUT_DIR" ]                                                    # -d teste l'existence du directory "<...>"
  then
    mkdir "$OUTPUT_DIR"
    Claim "$OUTPUT_DIR/ directory was not here yet : created" "/end"
  else
    Claim "$OUTPUT_DIR directory already created : ok" "/end"
  fi

  Build_lib
  Build_app
                                                                               # N.B. here the emptiness of $app_name is not tested because we never execute al the app : app_name is supposed to always be non null
  ./${app_name}.exe "${data_file%.nml}" < $DATA_DIR/$data_file > "$OUTPUT_DIR/"P1D_${app_name}_${data_file%.nml}_results.log # ${var%suffx} removes the suffix from the end of the var. "${data_file%.nml}" IS PASSED AS AN ARGUMENT OF THE APP ! /!\

  Claim "Done application"
  echo "################################################################################"
  echo "################################################################################"
  exit 0;;

*)
  Claim "something is weird : the command name is not recognized"
  echo "################################################################################"
  echo "################################################################################"
  exit 1;;
esac
