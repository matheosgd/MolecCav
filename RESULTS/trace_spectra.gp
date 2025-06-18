reset
set term qt font "Times, 12"                                                           
unset grid
show grid

#--- Full reso
File1 = '/home/segaud/MolecCav/OUT/Spectrum_1p1D_wmat1.88565999E-02_wcav1.88565999E-02_lamb3.00000003E-03.txt'
File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat1.88565999E-02_DTpm0._wcav1.88565999E-02_lamb3.00000003E-03.txt'
File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat1.88565999E-02_DTpm0._wcav1.88565999E-02_lamb3.00000003E-03.txt'
File4 = '/home/segaud/MolecCav/OUT/Spectrum_4p1D_wmat1.88565999E-02_wcav1.88565999E-02_lamb3.00000003E-03.txt'
#--- Cavity detuned
#File1 = '/home/segaud/MolecCav/OUT/Spectrum_1p1D_wmat1.88565999E-02_wcav1.89066008E-02_lamb3.00000003E-03.txt'
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat1.88565999E-02_DTpm0._wcav1.89066008E-02_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat1.88565999E-02_DTpm0._wcav1.89066008E-02_lamb3.00000003E-03.txt'
#--- Matter detuned DTmat = 2E-5
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat1.88565999E-02_DTpm1.99999995E-05_wcav1.88565999E-02_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat1.88565999E-02_DTpm1.99999995E-05_wcav1.88565999E-02_lamb3.00000003E-03.txt'
#--- Matter detuned DTmat = 5E-5
#File2 = '/home/segaud/MolecCav/OUT/TransInts_2p1D_wmat1.88565999E-02_DTpm1.99999995E-05_wcav1.88565999E-02_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat1.88565999E-02_DTpm4.99999987E-05_wcav1.88565999E-02_lamb3.00000003E-03.txt'
#--- Anar
File1_anar = '/home/segaud/MolecCav/OUT/Spectrum_ANAR_1p1D_wmat1.80664174E-02_wcav1.80664174E-02_lamb3.00000003E-03.txt'
File2_anar = '/home/segaud/MolecCav/OUT/Spectrum_ANAR_2p1D_wmat1.80664174E-02_DTpm0._wcav1.80664174E-02_lamb3.00000003E-03.txt'
#--- Not a spectrum
File = '/home/segaud/MolecCav/RESULTS/split.txt'

offset = 0 #coupling_strength*1E12

unset title #set title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'
set xlabel 'Transition energy [cm-1]'
set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]

unset xrange
unset yrange
#min = 5.8154906356422710E-003 - 1E-4
#max = 5.9170696435614355E-003 + 1E-4
#min = 4100
#max = 4180
#set xrange [min:max]
#set xrange [min*Conv:max*Conv]
#set yrange [0.0249:0.049]

set key left top

#  plot File1 u 1:2 w l lw 2 title '1p1D'
#replot File1_anar u 1:2 w l lw 2 title '1p1D (Anharmonic)'
#replot File2_anar u 1:2 w l lw 2 title '2p1D (Anharmonic)'
#replot File2 u 1:2 w l lw 2 title '2p1D'
#replot File3 u 1:2 w l lw 2 title '3p1D'
#replot File4 u 1:2 w l lw 2 title '4p1D'

set xlabel 'sqrt(N_{mat})'
set ylabel 'Band splitting [cm^{-1}]' #Coupling strength [a.u.]
Conv = 219474.6
stats '/home/segaud/MolecCav/RESULTS/split.txt' using (sqrt($1)):($2)*Conv #Cree les valeurs stats (STATS_correlation)
f(x) = a*x + b
set fit logfile "/home/segaud/MolecCav/RESULTS/split_stat.txt" #Indique dans quel fichier ecrire les resultats du fit
fit f(x) '/home/segaud/MolecCav/RESULTS/split.txt' u (sqrt($1)):($2)*Conv via a,b #le tracé
  plot File u (sqrt($1)):($2)*Conv w points pointtype 7 lw 2 title 'Splitting', f(x) t sprintf("ax+b R^2 = %g", STATS_correlation)

set key left top