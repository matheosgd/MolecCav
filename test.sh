#! /bin/bash

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Syntaxe :
# ./test.sh <name of the test (1)> <name of the data file (1) (optional)>
# (1) : without the extension
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

if [ -z "$1" ]                                             # z=!n; teste si la valeur est une chaîne vide.
then
  echo "##########################"
  echo "# No test name provided !#"
  echo "###### cf. Syntaxe #######"
  echo "##########################"
  exit 0
else 
  #echo "$1"
  test_name=$1
  echo "Running $test_name.f90 ..."
fi

if [ -n "$2" ]                                             # n teste si la valeur n'est pas une chaîne vide
then
  #echo "$2"
  data_file=$2
  echo "...using data of $data_file.nml"
else
  data_file="data_tests"
  echo "...using data of the DEFAULT data_tests.nml file"
fi

if [ ! -d "OBJ/" ]                                         # d teste l'existence du directory "<...>"
then
  mkdir OBJ/
  echo "OBJ/ directory was not here yet : created"
else
  echo "OBJ/ directory already created : ok"
fi

cd ~/MolecCav/OBJ/                                         # /!\This path may have to be changed /!\

gfortran -c ../SRC/MC_cavity_mode_m.f90
gfortran -c ../SRC/MC_operator_1D_m.f90
gfortran -c ../TESTS/${test_name}.f90
cd ~/MolecCav/                                             #/!\This path may have to be changed /!\
gfortran -o ${test_name}.exe Ext_Lib/QDUtilLib-main/libQD_gfortran_opt0_omp1_lapack1_int4.a OBJ/MC_cavity_mode_m.o OBJ/MC_operator_1D_m.o OBJ/${test_name}.o
./$test_name.exe < $data_file.nml > $data_file.out
