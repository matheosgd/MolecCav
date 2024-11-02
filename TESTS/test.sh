#! /usr/bin/bash

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Syntaxe :
# ./test.sh <name of the test without the extension>
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#echo "$1"
test_name=$1
#echo $test_name

cd ~/MolecCav/OBJ/ #/!\This path may have to be changed /!\

gfortran -c ../TESTS/${test_name}.f90
gfortran -o ${test_name}.exe ${test_name}.o

./$test_name.exe
