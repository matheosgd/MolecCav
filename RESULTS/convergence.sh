#! /bin/bash

rm -f "/home/msegaud/MolecCav/RESULTS/Eigenvalues_coupled_ham.txt"

cd ~/MolecCav
make all

coupling_strength="0.1"
for Nb in 3 4 5 6 7 8 9 10 11 12 #13 14 15 16 17 18
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
  lambda = $coupling_strength
  eq_pos = 0
  /
**
  echo "Finished Nb = $Nb"

  eigenvalues="$(grep "VP_TotH" OUT/App_MolecCav.log)"

  mid="$((1+${#eigenvalues}/2))"
  first_states="${eigenvalues:$mid:100}..."
  echo "$Nb ${first_states}" >> "/home/msegaud/MolecCav/RESULTS/Eigenvalues_coupled_ham.txt"

  echo "Done for Nb = $Nb"
done

cd /home/msegaud/MolecCav/RESULTS
gnuplot <<EOF
  load 'trace_cv.gp'
EOF
