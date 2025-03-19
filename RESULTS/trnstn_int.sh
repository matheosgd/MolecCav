#! /bin/bash

rm -f "/home/msegaud/MolecCav/RESULTS/transition_intensities.txt"

cd ~/MolecCav
make all
echo -e "\lambda -------- 0TO1 ------------ 0TO2 ------------ 0TO3 ------------ 0TO4" > "/home/msegaud/MolecCav/RESULTS/transition_intensities.txt"

for coupling_strength in 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24
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
  echo "Finished Nb = $coupling_strength"

  Intensities="$(grep "Intensities matrix   1" OUT/App_MolecCav.log)"
  GSto1="${Intensities:51: 12}"
  GSto2="${Intensities:70: 12}"
  GSto3="${Intensities:89: 12}"
  GSto4="${Intensities:108:12}"

  echo "$coupling_strength      ${GSto1}      ${GSto2}      ${GSto3}      ${GSto4}" >> "/home/msegaud/MolecCav/RESULTS/transition_intensities.txt"

  echo "Done for Nb = $coupling_strength"
done

#cd /home/msegaud/MolecCav/RESULTS
#gnuplot <<EOF                                   #pls trucs à tracer donc term png pas content donc on lance le script à la main
#  load 'trace_coup_eff.gp'
#EOF

echo -e "\nDetuning -------- 0TO1 ------------ 0TO2 ------------ 0TO3 ------------ 0TO4" >> "/home/msegaud/MolecCav/RESULTS/transition_intensities.txt"

for I in -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9
do 
  Detuning=$(echo "$I*0.0002" |bc -l)                                            #/!\ MANDATORY syntax to compute with floats !!!!! /!\ 
  w_C=$(echo "0.0058665+$Detuning" |bc -l)
  echo -e "\n Doing Detuning = $Detuning (w_C = $w_C)..."

  ./App_MolecCav.exe << ** > "OUT/App_MolecCav.log"
  &HO_1             !The diatomic molecule
  D = 1             !Label of the basis/HO/mode/dimension
  Nb = 10           !Number of basis vectors associated with the HO D
  w = 0.0058665     !Eigenpulsation associated with the HO D (HF molecule : 0.005866505831680149)
  m = 1744.60504565 !Mass associated with the HO D
  lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
  eq_pos = 0        !Equilibrium position of the HO
  /

  &HO_1             !The Cavity
  D = 2
  Nb = 11
  w = $w_C
  m = 1.0
  lambda = 0.1
  eq_pos = 0
  /
**
  echo "Finished Detuning = $Detuning"

  Intensities="$(grep "Intensities matrix   1" OUT/App_MolecCav.log)"
  GSto1="${Intensities:51: 12}"
  GSto2="${Intensities:70: 12}"
  GSto3="${Intensities:89: 12}"
  GSto4="${Intensities:108:12}"

  echo "$Detuning      ${GSto1}      ${GSto2}      ${GSto3}      ${GSto4}" >> "/home/msegaud/MolecCav/RESULTS/transition_intensities.txt"

  echo "Done for Detuning = ${Detuning}."
done

#cd /home/msegaud/MolecCav/RESULTS
#gnuplot <<EOF                                   #pls trucs à tracer donc term png pas content donc on lance le script à la main
#  load 'trace_coup_eff.gp'
#EOF
