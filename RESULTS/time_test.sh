#! /bin/bash

rm -f "/home/segaud/MolecCav/RESULTS/time/time_test_out.txt"
rm -f "/home/segaud/MolecCav/RESULTS/time/time_data_out.txt"

cd ~/MolecCav
make all MAIN=App_time

echo -e "Nb --------- time (cpu) [opt matrices] --------- time (cpu) [dense matrices]" > "/home/segaud/MolecCav/RESULTS/time/time_test_out.txt"

for Nb in 10 20 30 40 50
do 
  echo -e "\n Doing Nb = $Nb..."

  ./App_time.exe << ** > "OUT/App_time.log"
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
lambda = 0.04
eq_pos = 0
/

&HO_1             !The diatomic molecule
D = 3             !Label of the basis/HO/mode/dimension
Nb = $Nb          !Number of basis vectors associated with the HO D
w = 0.0058665     !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
m = 1744.60504565 !Mass associated with the HO D
lambda = 1.0      !Strength parameter of the coupling between the mode D and the molecule
eq_pos = 0        !Equilibrium position of the HO
/

&HO_1             !The Cavity
D = 4
Nb = $Nb
w = 0.0068665
m = 1.0
lambda = 0.04
eq_pos = 0
/
**
  echo "Finished computation for Nb = $Nb"

  time="$(grep "Total cpu (s):" OUT/App_time.log)"
  echo -e "----------Nb = $Nb----------\n$time\n"                               >> "/home/segaud/MolecCav/RESULTS/time/time_data_out.txt"

  time_opt="${time:62:9}"
  time_dense="$(echo "${time:98:9} - $time_opt" |bc -l)"
  echo "$Nb           ${time_opt}                                ${time_dense}" >> "/home/segaud/MolecCav/RESULTS/time/time_test_out.txt"

  echo "Done for Nb = $Nb"
done

cd /home/segaud/MolecCav/RESULTS/time

gnuplot -p -e "end_nb=50; end_pt=5" /home/segaud/MolecCav/RESULTS/time/trace_time.gp
gnuplot -p -e "end_nb=40; end_pt=4" /home/segaud/MolecCav/RESULTS/time/trace_time.gp