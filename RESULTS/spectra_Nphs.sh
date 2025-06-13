#! /bin/bash

Nb="10"
MAIN="Spec_Nphs" #App_trnstn_int App_trnstn_int_ExactSmllBss

if [ ! -d "/home/segaud/MolecCav/RESULTS/spectra/" ]                                               # d teste l'existence du directory "<...>"
then
  mkdir "/home/segaud/MolecCav/RESULTS/spectra/"
else
  rm -f "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
fi

cd ~/MolecCav
make all MAIN="$MAIN"
echo -e "reset"                                                                                     > "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "set term qt font \"Times, 12\""                                                           >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "unset grid\nshow grid"                                                                    >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "\nunset title"                                                                              >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "set xlabel 'Transition energy [cm-1]'"                                                    >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]"            >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "\nunset xrange\nunset yrange"                                                             >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "#min = 5.7455455693706022E-003 - 1E-4"                                                    >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "#max = 5.9850104937492297E-003 + 1E-4"                                                    >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "#set xrange [min*Conv:max*Conv]"                                                          >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "#set yrange [0.0249:0.049]"                                                               >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
echo -e "\n#set key left top"                                                                       >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"

for N_phs in 1 2 3
do 
  echo -e "\n Doing number of photons = $N_phs..."
  w_mat=$(echo "0.0058665*$N_phs" |bc -l)
  echo -e "... i.e. doing matter mode frequency = $w_mat..."

  ./${MAIN}.exe << ** > "OUT/${MAIN}${N_phs}.log"
  &TOTAL_HAMILTONIAN
  N_mat = 1
  N_cav = 1
  /

  !--- NEEED ONLY ONE NML FOR ALL OPERATORS (related to the same modes) !!!
  &Matter_mode
  lambda =  0.05477225575051661 ! 6.32e-2   !Strength parameter of the coupling between the mode D and the molecule
  CoeffsDipMomt = 0.0, 1.0, 0.0,  0.0,  0.0 !The n-derivatives of the matter dipole moment with respect to this mode's DOF i.e. the coefficients of its Taylor expansion
  Nb = 5                                    !Number of basis vectors associated with the HO D
  w = ${w_mat}                              !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
  m = 1744.60504565                         !Mass associated with the HO D
  Nq = 1
  Eq_pos = 0.0                              !Equilibrium position of the HO
  Scale_q = 0.0
  /

  &Cavity_mode
  lambda =  0.05477225575051661 ! 6.32e-2   !Strength parameter of the coupling between the mode D and the molecule
  Nb = 5                                    !Number of basis vectors associated with the HO D
  w = 0.0058665                             !Eigenpulsation associated with the HO D (maybe HF molecule : 0.005866505831680149)
  m = 1.0                                   !Mass associated with the HO D
  Nq = 1
  Eq_pos = 0.0                              !Equilibrium position of the HO
  Scale_q = 0.0
  /

  &DIPOLE_MOMENT !Please keep the same parameters as in total_hamiltonian, this point will be simplified later
  N_mat = 1
  N_cav = 1
  /
**
  echo "Finished computation for number of photons = $N_phs."

  File="$(grep "OUT/Spectrum_Nphs_" OUT/${MAIN}${N_phs}.log)"
  cd ~/MolecCav
  # echo "${PWD}/${File:1}"
  echo -e "\nFile${N_phs} = '${PWD}/${File:1}'"                                                    >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
  echo -e "offset${N_phs} = 0"                                                                     >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"  

  if [ "$N_phs" == 1 ]; then 
    echo -e "  plot File${N_phs} u 1:2 w l lw 2 title '${N_phs} photon(s)'"                        >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
    # echo -e "  plot File${N_phs} u 1:2 w l lw 2 notitle"                                           >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
  elif [ "$N_phs" == 2 ]; then 
    echo -e "replot File${N_phs} u 1:2 w l lw 2 title '${N_phs} photon(s)'"                        >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
    # echo -e "replot File${N_phs} u 1:2 w l lw 2 notitle"                                           >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
  elif [ "$N_phs" == 3 ]; then 
    echo -e "replot File${N_phs} u 1:2 w l lw 2 title '${N_phs} photon(s)'"                        >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
    # echo -e "replot File${N_phs} u 1:2 w l lw 2 notitle"                                           >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
  else 
    echo -e "replot File${N_phs} u 1:2 w l lw 2 title '${N_phs} photon(s)'"                        >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"
  fi

  echo -e "\n Done testing number of photons"
done

echo -e "\n#set key left top"                                                                      >> "/home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp"

gnuplot -s -p /home/segaud/MolecCav/RESULTS/spectra/trace_spectra_Nphs.gp