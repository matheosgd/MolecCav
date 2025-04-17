#! /bin/bash

rm -f "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_detu_out.txt"

cd ~/MolecCav
make all MAIN=App_perturb_effect
make all MAIN=Calc_perturb
#echo "DT ------ E1 ------------------------------- E2 ------------------------------- A(DT) ---------------------------- DT/A(0)" > "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_detu_out.txt"
echo "DT ------ E1 ------------------------------- E2 ---------------------------- E1_pert --------------- E2_pert --------------- DT/A(0)" > "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_detu_out.txt"

coupling_strength="0.04"
for DT in 0.00000 0.00002 0.00004 0.00006 0.00008 0.00010 0.00012 0.00014 0.00016 0.00018 0.00020 0.00022 0.00024 0.0002 0.0004 0.0006 0.0008 0.0010
do 
  echo -e "\n Doing detuning = $DT..."
  w_C=$(echo "0.0058665+$DT" |bc -l)

  ./App_perturb_effect.exe << ** > "OUT/App_perturb_effect_in_detu.log"
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
  w = ${w_C}
  m = 1.0
  lambda = ${coupling_strength}
  eq_pos = 0
  /
**
  echo "Finished computation for detuning = $DT"

  Ratio_all="$(grep "xxx DT" OUT/App_perturb_effect_in_detu.log)"
  Ratio="${Ratio_all:34}"
  Energy_level_1="$(grep "Energy level of the first excited state" OUT/App_perturb_effect_in_detu.log)"
  E1="${Energy_level_1:43}"
  Energy_level_2="$(grep "Energy level of the secnd excited state" OUT/App_perturb_effect_in_detu.log)"
  E2="${Energy_level_2:43}"
#  Gap_all="$(grep "2xA(DT)" OUT/App_perturb_effect_in_detu.log)"
#  Gap="${Gap_all:14}"

  ./Calc_perturb.exe << ** > "OUT/Calc_perturb_effect_detu.log"
  &PARAMETERS
  DT_nml = $DT
  lambda_nml = 0.04
  Case_nml = 1         !Means regime A >> DT
  /
**
  echo "Finished perturbation for detuning = $DT"

  Energy_levels="$(grep "Energy levels RC" OUT/Calc_perturb_effect_detu.log)"
  E1_pert="${Energy_levels:49:13}"
  E2_pert="${Energy_levels:68}"

#  echo "${DT}   ${E1}           ${E2}           ${Gap}           ${Ratio}" >> "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_detu_out.txt"
  echo "${DT}   ${E1}           ${E2}       ${E1_pert}           ${E2_pert}           ${Ratio}" >> "/home/segaud/MolecCav/RESULTS/perturbations/perturb_effect_in_detu_out.txt"

  echo "Done for detuning = $DT"
done

gnuplot -p -e "lambda=${coupling_strength}" /home/segaud/MolecCav/RESULTS/perturbations/trace_perturb_effect_in_detu.gp
