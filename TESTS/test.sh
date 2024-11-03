#! /usr/bin/bash

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

cd ~/MolecCav/OBJ/                                         #/!\This path may have to be changed /!\

if [ -n "$2" ]
then
  #echo "$2"
  data_file=$2
  echo "...using data of $data_file.nml"
else
  data_file="data"
  echo "...using data of the DEFAULT data.nml file"
fi

gfortran -c ../TESTS/${test_name}.f90
cd ~/MolecCav/                                             #/!\This path may have to be changed /!\
gfortran -o ${test_name}.exe OBJ/${test_name}.o
./$test_name.exe < $data_file.nml
