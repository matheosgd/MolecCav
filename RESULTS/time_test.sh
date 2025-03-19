#! /bin/bash

rm -f "/home/msegaud/MolecCav/RESULTS/time_test.txt"

cd ~/MolecCav
make all
echo -e "\n Nb --------- time (cpu) [dense matrices]" >> "/home/msegaud/MolecCav/RESULTS/time_test.txt"

for Nb in 10 20 30 40 50
do 
  echo -e "\n Doing Nb = $Nb..."

  ./App_MolecCav.exe << ** > "OUT/App_MolecCav.log"
&HO_1             !The diatomic molecule
D = 1             !Label of the basis/HO/mode/dimension
Nb = $Nb          !Number of basis vectors associated with the HO D
w = 0.0058665     !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
m = 1744.60504565 !Mass associated with the HO D
lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
eq_pos = 0        !Equilibrium position of the HO
/

&HO_1             !The Cavity
D = 2
Nb = $Nb
w = 0.0068665
m = 1.0
lambda = 0.1
eq_pos = 0
/
**
  echo "Finished Nb = $Nb"

  time="$(grep "Total cpu (s):" OUT/App_MolecCav.log)"
  time_cut="${time:65} ..."

  echo "$Nb ${time_cut}" >> "/home/msegaud/MolecCav/RESULTS/time_test.txt"

  echo "Done for Nb = $Nb"
done

cd /home/msegaud/MolecCav/RESULTS
#gnuplot <<EOF
#  load 'trace_time.gp'
#EOF
