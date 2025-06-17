reset
set term qt font "Times, 10"                                                           
unset grid
show grid

#--- Full reso
File1 = '/home/segaud/MolecCav/OUT/Spectrum_1p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_DTpm0._wcav5.86650008E-03_lamb3.00000003E-03.txt'
File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_DTpm0._wcav5.86650008E-03_lamb3.00000003E-03.txt'
File4 = '/home/segaud/MolecCav/OUT/Spectrum_4p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#--- Cavity detuned
#File1 = '/home/segaud/MolecCav/OUT/Spectrum_1p1D_wmat5.86650008E-03_wcav6.06650021E-03_lamb3.00000003E-03.txt'
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_DTpm0._wcav6.06650021E-03_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_DTpm0._wcav6.06650021E-03_lamb3.00000003E-03.txt'
#--- Matter detuned DTmat = 1E-4
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_DTpm9.99999975E-05_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_DTpm9.99999975E-05_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#--- Matter detuned DTmat = 2E-4
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_DTpm1.99999995E-04_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_DTpm1.99999995E-04_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#--- Anar
#File1 = '/home/segaud/MolecCav/OUT/Spectrum_ANAR_1p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_ANAR_2p1D_wmat5.86650008E-03_DTpm9.99999975E-05_wcav5.86650008E-03_lamb3.00000003E-03.txt'


offset = 0 #coupling_strength*1E12

unset title #set title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'

unset xrange
unset yrange
#min = 5.8154906356422710E-003 - 1E-4
#max = 5.9170696435614355E-003 + 1E-4
#set xrange [min:max]
#set xrange [min*Conv:max*Conv]
#set yrange [0.0249:0.049]

#set key left top

set yrange [0:2.1]
set zrange [0:0.035]
 splot File1 u 1:3:2 w l lw 2 title '1p1D'
replot File2 u 1:3:2 w l lw 2 title '2p1D'
replot File3 u 1:3:2 w l lw 2 title '3p1D'


#set multiplot
#set size 1.0,0.333
#set origin 0.0,0.666
#unset xlabel
#set ylabel ' '
#  plot File1 u 1:2 w l lw 2 title '1p1D'
#set size 1.0,0.333
#set origin 0.0,0.333
#set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]
#  plot File2 u 1:2 w l lw 2 title '2p1D'
#set size 1.0,0.333
#set origin 0.0,0.0
#set ylabel ' '
#set xlabel 'Transition energy [cm-1]'
#  plot File3 u 1:2 w l lw 2 title '3p1D'
#set size 0.4,0.4
#set origin 0.1,0.1
#plot File4 u 1:2 w l lw 2 title '4p1D'
#unset multiplot


#set key left top