#! /bin/bash

rm -f "/home/segaud/MolecCav/RESULTS/perturbations.perturbations_sumup.txt"

cd ~/MolecCav
make all MAIN=App_perturbations
make all MAIN=Calc_perturb


./App_perturbations.exe << ** > "OUT/App_perturbations_smallB.log"
&HO_1             !The diatomic molecule
D = 1             !Label of the basis/HO/mode/dimension
Nb = 2            !Number of basis vectors associated with the HO D
w = 0.0058        !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
m = 1744.60504565 !Mass associated with the HO D
lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
eq_pos = 0        !Equilibrium position of the HO
/

&HO_1             !The Cavity_RuC : \lambda = DT = 0 = Resonant/unCoupled (a bit useless actually => commented hereafter)
D = 2
Nb = 2
w = 0.0058
m = 1.0
lambda = 0.00
eq_pos = 0
/

&HO_1             !The Cavity_RC : A(DT) >> DT : zeroth order = Resonant/Coupled
D = 3
Nb = 2
w = 0.00582       !Perturbation = Detuning, here DT = 2e-5
m = 1.0
lambda = 0.04     !Zeroth order = Resonant/Coupled
eq_pos = 0
/

&HO_1             !The Cavity_oRnC : DT >> A(DT) : zeroth order = off-Resonance/not Coupled
D = 4
Nb = 2
w = 0.0060        !Zeroth order = off-Resonance/not Coupled, here DT = 2e-4
m = 1.0
lambda = 0.004    !Perturbation = Coupling
eq_pos = 0
/
**
echo "Finished App_perturbations_smallB"

./App_perturbations.exe << ** > "OUT/App_perturbations_convergenceB.log"
&HO_1             !The diatomic molecule
D = 1             !Label of the basis/HO/mode/dimension
Nb = 10            !Number of basis vectors associated with the HO D
w = 0.0058        !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
m = 1744.60504565 !Mass associated with the HO D
lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
eq_pos = 0        !Equilibrium position of the HO
/

&HO_1             !The Cavity_RuC : \lambda = DT = 0 = Resonant/unCoupled
D = 2
Nb = 11
w = 0.0058
m = 1.0
lambda = 0.00
eq_pos = 0
/

&HO_1             !The Cavity_RC : A(DT) >> DT : zeroth order = Resonant/Coupled
D = 2
Nb = 11
w = 0.00582       !Perturbation = Detuning, here DT = 2e-5
m = 1.0
lambda = 0.008    !Zeroth order = Resonant/Coupled
eq_pos = 0
/

&HO_1             !The Cavity_oRnC : DT >> A(DT) : zeroth order = off-Resonance/not Coupled
D = 3
Nb = 11
w = 0.0060        !Zeroth order = off-Resonance/not Coupled, here DT = 2e-4
m = 1.0
lambda = 0.004    !Perturbation = Coupling
eq_pos = 0
/
**
echo "Finished App_perturbations_convergenceB"

./Calc_perturb.exe << ** > "OUT/Calc_perturb.log"
&PARAMETERS
DT_nml = $DT
lambda_nml = 0.04
Case_nml = 0         !Means the nml will not be used
/
**
echo "Finished Calc_perturb"

cp "OUT/Calc_perturb.log" "RESULTS/perturbations/perturbations_sumup.txt"

echo -e "\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxSmall Basis setxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"        >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"        >> "RESULTS/perturbations/perturbations_sumup.txt"
#Energy_levels="$(grep "Energy levels TotH_RuC" OUT/App_perturbations_smallB.log)"
#Gap="$(grep "Energy gap : |E_2 - E_1| (RuC)" OUT/App_perturbations_smallB.log)"
#echo "${Energy_levels::81}" >> "RESULTS/perturbations/perturbations_sumup.txt"
#echo "$Gap"                 >> "RESULTS/perturbations/perturbations_sumup.txt"

Energy_levels="$(grep "Energy levels TotH_RC" OUT/App_perturbations_smallB.log)"
Gap="$(grep "Energy gap : |E_2 - E_1| (RC)" OUT/App_perturbations_smallB.log)"
echo -e "\n${Energy_levels::81}" >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "$Gap"                      >> "RESULTS/perturbations/perturbations_sumup.txt"

Energy_levels="$(grep "Energy levels TotH_oRuC" OUT/App_perturbations_smallB.log)"
Gap="$(grep "Energy gap : |E_2 - E_1| (oRuC)" OUT/App_perturbations_smallB.log)"
echo -e "\n${Energy_levels::81}" >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "$Gap"                      >> "RESULTS/perturbations/perturbations_sumup.txt"

echo -e "\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxConverged Basis setxxxxxxxxxxxxxxxxxxxxxxxxxxxx"        >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"        >> "RESULTS/perturbations/perturbations_sumup.txt"
#Energy_levels="$(grep "Energy levels TotH_RuC" OUT/App_perturbations_convergenceB.log)"
#Gap="$(grep "Energy gap : |E_2 - E_1| (RuC)" OUT/App_perturbations_convergenceB.log)"
#echo "${Energy_levels::81}" >> "RESULTS/perturbations/perturbations_sumup.txt"
#echo "$Gap" >> "RESULTS/perturbations/perturbations_sumup.txt"

Energy_levels="$(grep "Energy levels TotH_RC" OUT/App_perturbations_convergenceB.log)"
Gap="$(grep "Energy gap : |E_2 - E_1| (RC)" OUT/App_perturbations_convergenceB.log)"
echo -e "\n${Energy_levels::81}" >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "$Gap" >> "RESULTS/perturbations/perturbations_sumup.txt"

Energy_levels="$(grep "Energy levels TotH_oRuC" OUT/App_perturbations_convergenceB.log)"
Gap="$(grep "Energy gap : |E_2 - E_1| (oRuC)" OUT/App_perturbations_convergenceB.log)"
echo -e "\n${Energy_levels::81}" >> "RESULTS/perturbations/perturbations_sumup.txt"
echo "$Gap" >> "RESULTS/perturbations/perturbations_sumup.txt"

#execute the perturb_effect_*.sh