reset
set term qt font "Times, 12"                                                           
unset grid
show grid

#--- Full reso
File1 = '/home/segaud/MolecCav/OUT/Spectrum_1p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
#--- Cavity detuned
#File1 = '/home/segaud/MolecCav/OUT/Spectrum_1p1D_wmat5.86650008E-03_wcav6.06650021E-03_lamb3.00000003E-03.txt'
#File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_wcav6.06650021E-03_lamb3.00000003E-03.txt'
#File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_wcav6.06650021E-03_lamb3.00000003E-03.txt'
#--- Matter detuned
File2 = '/home/segaud/MolecCav/OUT/Spectrum_2p1D_wmat5.86650008E-03_DTpm9.99999975E-05_wcav5.86650008E-03_lamb3.00000003E-03.txt'
File3 = '/home/segaud/MolecCav/OUT/Spectrum_3p1D_wmat5.86650008E-03_DTpm9.99999975E-05_wcav5.86650008E-03_lamb3.00000003E-03.txt'

offset = 0 #coupling_strength*1E12

unset title #set title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'
set xlabel 'Transition energy [Ha]'
set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]

unset xrange
unset yrange
#min = 5.8154906356422710E-003 - 1E-4
#max = 5.9170696435614355E-003 + 1E-4
#set xrange [min:max]
#set xrange [min*Conv:max*Conv]
#set yrange [0.0249:0.049]

#set key left top

  plot File1 u 1:2 w l lw 2 title '1p1D'
replot File2 u 1:2 w l lw 2 title '2p1D'
replot File3 u 1:2 w l lw 2 title '3p1D'

#set key left top