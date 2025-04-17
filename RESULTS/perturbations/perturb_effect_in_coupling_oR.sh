#! /bin/bash

rm -f "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_coupling_oR_out.txt"

cd ~/MolecCav
make all MAIN=App_perturb_effect
make all MAIN=Calc_perturb
#echo "lambda -- E1 ------------------------------- E2 ------------------------------- A(DT) ---------------------------- DT/A(0)" > "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_coupling_oR_out.txt"
echo "lambda -- E1 ------------------------------- E2 ---------------------------- E1_pert --------------- E2_pert --------------- DT/A(DT)" > "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_coupling_oR_out.txt"

DT="0.0002" #DT << lambda but A(DT) >> lambda, A(DT)  being the relevant parameter to compare which one pertubates the other (cf matrices in Calc_perturb.f90/log)
w_C=$(echo "0.0058665+$DT" |bc -l)
for coupling_strength in 0.000 0.0015 0.003 0.0045 0.006 0.0075 0.0090 0.0105 0.0120 0.0135 0.0150 0.0165 0.0180 0.0195 0.02 #0.022 0.024 0.02 0.04 0.05
do 
  echo -e "\n Doing coupling strength = $coupling_strength..."

  ./App_perturb_effect.exe << ** > "OUT/App_perturb_effect_in_coupling_oR.log"
  &HO_1             !The diatomic molecule
  D = 1             !Label of the basis/HO/mode/dimension
  Nb = 2            !Number of basis vectors associated with the HO D
  w = 0.0058665     !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
  m = 1744.60504565 !Mass associated with the HO D
  lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
  eq_pos = 0        !Equilibrium position of the HO
  /

  &HO_1             !The Cavity
  D = 2
  Nb = 2
  w = $w_C
  m = 1.0
  lambda = $coupling_strength
  eq_pos = 0
  /
**
  echo "Finished computation for coupling strength = $coupling_strength"

  Ratio_all="$(grep "xxx DT" OUT/App_perturb_effect_in_coupling_oR.log)"
  Ratio="${Ratio_all:34}"
  Energy_level_1="$(grep "Energy level of the first excited state" OUT/App_perturb_effect_in_coupling_oR.log)"
  E1="${Energy_level_1:43}"
  Energy_level_2="$(grep "Energy level of the secnd excited state" OUT/App_perturb_effect_in_coupling_oR.log)"
  E2="${Energy_level_2:43}"
#  Gap_all="$(grep "2xA(DT)" OUT/App_perturb_effect_in_coupling_oR.log)"
#  Gap="${Gap_all:14}"

  ./Calc_perturb.exe << ** > "OUT/Calc_perturb_effect_coupling_oR.log"
  &PARAMETERS
  DT_nml = $DT
  lambda_nml = $coupling_strength
  Case_nml = 2         !Means regime DT >> A
  /
**
  echo "Finished perturbation for coupling strength = $coupling_strength"

  Energy_levels="$(grep "Energy levels hRnC" OUT/Calc_perturb_effect_coupling_oR.log)"
  E1_pert="${Energy_levels:51:13}"
  E2_pert="${Energy_levels:70}"

#  echo "${DT}   ${E1}           ${E2}           ${Gap}           ${Ratio}" >> "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_coupling_oR_out.txt"
  echo "${coupling_strength}   ${E1}           ${E2}       ${E1_pert}           ${E2_pert}           ${Ratio}" >> "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_coupling_oR_out.txt"

  echo "Done for coupling strength = $coupling_strength"
done

gnuplot -p -e "DT=${DT}" /home/segaud/MolecCav/RESULTS/perturbations/trace_perturb_effect_in_coupling_oR.gp
