#! /bin/bash

rm -f "/home/msegaud/MolecCav/RESULTS/coupling_effect.txt"

cd ~/MolecCav
make all
echo -e "\n lambda --------- w1 ---------------- w2" >> "/home/msegaud/MolecCav/RESULTS/coupling_effect.txt"

for coupling_strength in 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24
do 
  echo -e "\n Doing Nb = $coupling_strength..."

  ./App_MolecCav.exe << ** > "OUT/App_MolecCav.log"
  &HO_1             !The diatomic molecule
  D = 1             !Label of the basis/HO/mode/dimension
  Nb = 10           !Number of basis vectors associated with the HO D
  w = 0.0058665     !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
  m = 1744.60504565 !Mass associated with the HO D
  lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
  eq_pos = 0        !Equilibrium position of the HO
  /

  &HO_1             !The Cavity
  D = 2
  Nb = 11
  w = 0.0068665
  m = 1.0
  lambda = $coupling_strength
  eq_pos = 0
  /
**
  echo "Finished Nb = $coupling_strength"

  w1_all="$(grep "1^{th} Normal coordinate" OUT/App_MolecCav.log)"
  w2_all="$(grep "2^{th} Normal coordinate" OUT/App_MolecCav.log)"

  w1="${w1_all: -24}" #/!\ space needed here !
  w2="${w2_all: -24}"

  echo "$coupling_strength ${w1} ${w2}" >> "/home/msegaud/MolecCav/RESULTS/coupling_effect.txt"

  echo "Done for Nb = $coupling_strength"
done

cd /home/msegaud/MolecCav/RESULTS
#gnuplot <<EOF                                   #pls trucs à tracer donc term png pas content donc on lance le script à la main
#  load 'trace_coup_eff.gp'
#EOF
