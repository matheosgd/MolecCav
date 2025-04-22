#! /bin/bash

if [ ! -d "/home/segaud/MolecCav/RESULTS/int_strength/" ]                                                             # d teste l'existence du directory "<...>"
then
  mkdir "/home/segaud/MolecCav/RESULTS/int_strength/"
else
  rm -f "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
fi

cd ~/MolecCav
make all MAIN=App_trnstn_int
echo -e "reset"                                                                                    >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "set term qt font \"Times, 12\""                                                           >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "set grid\nshow grid"                                                                      >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "\nset title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'" >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "set xlabel 'Transition energy [Ha]'"                                                      >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "set ylabel 'Arbitrary units' #Coupling strength [a.u.]"                                   >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "\nunset xrange\nunset yrange"                                                             >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "min = 5.7455455693706022E-003 - 1E-4"                                                     >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "max = 5.9850104937492297E-003 + 1E-4"                                                     >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "set xrange [min:max]"                                                                     >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "#set yrange [0.0249:0.049]"                                                               >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "#set key left top"                                                                        >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "\nConv = 21947.46 # 1Ha = Conv.cm-1"                                                      >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "Gam = 0.1/Conv  # 10cm-1"                                                                 >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "set samples 800\nshow samples"                                                            >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
echo -e "L(x, x_0, Gam) = ( Gam/(2*pi) ) / ( ((Gam**2)/4) + (x-x_0)**2 )"                          >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"

for coupling_strength in 0.00 0.002 0.004 0.006 0.008 0.01
do 
  echo -e "\n Doing coupling_strength = $coupling_strength..."

  ./App_trnstn_int.exe << ** > "OUT/App_trnstn_int.log"
  &HO_1                       !The diatomic molecule
  D = 1                       !Label of the basis/HO/mode/dimension
  Nb = 10                     !Number of basis vectors associated with the HO D
  w = 0.0058665               !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
  m = 1744.60504565           !Mass associated with the HO D
  lambda = $coupling_strength !Strength parameter of the coupling between the mode D and the molecule
  eq_pos = 0                  !Equilibrium position of the HO
  /

  &HO_1                       !The Cavity
  D = 2
  Nb = 11
  w = 0.0058665
  m = 1.0
  lambda = $coupling_strength
  eq_pos = 0
  /
**
  echo "Finished computation for coupling strength = $coupling_strength"

  Intensities="$(grep "Intensities matrix   1" OUT/App_trnstn_int.log)"
  GSto1="${Intensities:51: 12}"
  GSto2="${Intensities:70: 12}"
  Enrgy1="$(grep "Transition energy GSto1" OUT/App_trnstn_int.log)"
  Enrgy1="${Enrgy1:27}"
  Enrgy2="$(grep "Transition energy GSto2" OUT/App_trnstn_int.log)"
  Enrgy2="${Enrgy2:27}"
#  GSto3="${Intensities:89: 12}"
#  GSto4="${Intensities:108:12}"
#  Enrgy3="$(grep "Transition energy GSto3" OUT/App_trnstn_int.log)"
#  Enrgy3="${Enrgy3:27}"
#  Enrgy4="$(grep "Transition energy GSto4" OUT/App_trnstn_int.log)"
#  Enrgy4="${Enrgy4:27}"

  echo -e "\noffset = 0 #$coupling_strength*1E12" >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"  
  if [ "$coupling_strength" == 0.00 ]; then 
    echo -e "  plot ( ${GSto1}*L(x, ${Enrgy1}, Gam) + ${GSto2}*L(x, ${Enrgy2}, Gam) + offset) w l lw 2 t '\lambda = 0'" >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
  else 
    echo "replot ( ${GSto1}*L(x, ${Enrgy1}, Gam) + ${GSto2}*L(x, ${Enrgy2}, Gam) + offset) w l lw 2 t '\lambda = $coupling_strength'" >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"
  fi

  echo -e "\n Done testing coupling strength"
done

echo -e "\n#set key left top" >> "/home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp"

gnuplot -p /home/segaud/MolecCav/RESULTS/int_strength/trace_int_strength_out.gp