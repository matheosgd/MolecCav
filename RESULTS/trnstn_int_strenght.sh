#! /bin/bash

rm -f "/home/msegaud/MolecCav/RESULTS/transition_intensities.txt"

cd ~/MolecCav
make all
echo -e "lambda --- Transition energy 0to1 --- Transition intensity 0to1 --- Transition energy 0to2 --- Transition intensity 0to2" > "/home/segaud/MolecCav/RESULTS/transition_intensities.txt"

for coupling_strength in 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10
do 
  echo -e "\n Doing coupling_strength = $coupling_strength..."

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
  w = 0.0058665
  m = 1.0
  lambda = $coupling_strength
  eq_pos = 0
  /
**
  echo "Finished coupling strength = $coupling_strength"

  Intensities="$(grep "Intensities matrix   1" OUT/App_MolecCav.log)"
  GSto1="${Intensities:51: 12}"
  GSto2="${Intensities:70: 12}"
#  GSto3="${Intensities:89: 12}"
#  GSto4="${Intensities:108:12}"
  Enrgy1="$(grep "Transition energy GSto1" OUT/App_MolecCav.log)"
  Enrgy1="${Enrgy1:27}"
  Enrgy2="$(grep "Transition energy GSto2" OUT/App_MolecCav.log)"
  Enrgy2="${Enrgy2:27}"
#  Enrgy3="$(grep "Transition energy GSto3" OUT/App_MolecCav.log)"
#  Enrgy3="${Enrgy3:27}"
#  Enrgy4="$(grep "Transition energy GSto4" OUT/App_MolecCav.log)"
#  Enrgy4="${Enrgy4:27}"

  echo "$coupling_strength       ${Enrgy1}    ${GSto1}                  ${Enrgy2}    ${GSto2}" >> "/home/segaud/MolecCav/RESULTS/transition_intensities.txt"

  echo "Done testing coupling strength"
done

#cd /home/msegaud/MolecCav/RESULTS
#gnuplot <<EOF                                   #pls trucs à tracer donc term png pas content donc on lance le script à la main
#  load 'trace_coup_eff.gp'
#EOF