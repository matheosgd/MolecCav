#! /bin/bash

if [ ! -d "/home/segaud/MolecCav/RESULTS/int_detu/" ]                                                             # d teste l'existence du directory "<...>"
then
  mkdir "/home/segaud/MolecCav/RESULTS/int_detu/"
else
  rm -f "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
fi

cd ~/MolecCav
make all MAIN=App_trnstn_int
echo -e "reset"                                                                                                 >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "set term qt font \"Times, 12\""                                                                        >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "set grid\nshow grid"                                                                                   >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "\nset title  'Intensities = f(detuning) for small coupling (\lambda = 4e-3 a.u.; w_M = w_{HF}) [a.u.]" >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "set xlabel 'Transition energy [Ha]'"                                                                   >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "set ylabel 'Arbitrary units [a.u.]'"                                                                   >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "\nunset xrange\nunset yrange"                                                                          >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "min = 5.4610394673397069E-003 - 1E-5"                                                                  >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "max = 6.2723303490619256E-003 + 1E-5"                                                                  >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "set xrange [min:max]"                                                                                  >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "#set yrange [0.0249:0.049]"                                                                            >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "#set key left top"                                                                                     >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "\nConv = 21947.46 # 1Ha = Conv.cm-1"                                                                   >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "Gam = 0.3/Conv  # 30cm-1"                                                                              >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "set samples 999\nshow samples"                                                                         >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
echo -e "L(x, x_0, Gam) = ( Gam/(2*pi) ) / ( ((Gam**2)/4) + (x-x_0)**2 )"                                       >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"

for I in 0 1 -1 #2 -2
do 
  Detuning=$(echo "$I*0.0002" |bc -l)                                          #/!\ MANDATORY syntax to compute with floats !!!!! /!\ 
  w_C=$(echo "0.0058665+$Detuning" |bc -l)
  echo -e "\n Doing Detuning = $Detuning (w_C = $w_C)..."

  ./App_trnstn_int.exe << ** > "OUT/App_trnstn_int.log"
  &HO_1             !The diatomic molecule
  D = 1             !Label of the basis/HO/mode/dimension
  Nb = 10           !Number of basis vectors associated with the HO D
  w = 0.0058665     !Eigenpulsation associated with the HO D (HF molecule : 0.005866505831680149)
  m = 1744.60504565 !Mass associated with the HO D
  lambda = 0.004    !Strength parameter of the coupling between the mode D and the molecule
  eq_pos = 0        !Equilibrium position of the HO
  /

  &HO_1             !The Cavity
  D = 2
  Nb = 11
  w = $w_C
  m = 1.0
  lambda = 0.004
  eq_pos = 0
  /
**
  echo "Finished computation for detuning = $Detuning"

  Intensities="$(grep "Intensities matrix   1" OUT/App_trnstn_int.log)"
  GSto1="${Intensities:51: 12}"
  GSto2="${Intensities:70: 12}"
  Enrgy1="$(grep "Transition energy GSto1" OUT/App_trnstn_int.log)"
  Enrgy1="${Enrgy1:27}"
  Enrgy2="$(grep "Transition energy GSto2" OUT/App_trnstn_int.log)"
  Enrgy2="${Enrgy2:27}"
#  GSto3="${Intensities:89: 12}"
#  GSto4="${Intensities:108:12}"

  echo -e "\noffset = 0 #$Detuning*1E12" >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"  
  if [ "$I" == "0" ]; then 
    echo -e "  plot ( ${GSto1}*L(x, ${Enrgy1}, Gam) + ${GSto2}*L(x, ${Enrgy2}, Gam) + offset) w l lw 2 t 'Detuning = $Detuning'" >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
  else 
    echo "replot ( ${GSto1}*L(x, ${Enrgy1}, Gam) + ${GSto2}*L(x, ${Enrgy2}, Gam) + offset) w l lw 2 t 'Detuning = $Detuning'"    >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"
  fi

  echo -e "\n Done testing coupling strength"
done

echo -e "#\nset key left top" >> "/home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp"

gnuplot -p /home/segaud/MolecCav/RESULTS/int_detu/trace_int_detu_out.gp