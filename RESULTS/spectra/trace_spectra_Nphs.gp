reset
set term qt font "Times, 12"
unset grid
show grid

unset title
set xlabel 'Transition energy [cm-1]'
set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]

unset xrange
unset yrange
#min = 5.7455455693706022E-003 - 1E-4
#max = 5.9850104937492297E-003 + 1E-4
#set xrange [min*Conv:max*Conv]
#set yrange [0.0249:0.049]

#set key left top

File1 = '/home/segaud/MolecCav/OUT/Spectrum_Nphs_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
offset1 = 0
  plot File1 u 1:2 w l lw 2 notitle

File2 = '/home/segaud/MolecCav/OUT/Spectrum_Nphs_wmat1.17330002E-02_wcav5.86650008E-03_lamb3.00000003E-03.txt'
offset2 = 0
replot File2 u 1:2 w l lw 2 notitle

File3 = '/home/segaud/MolecCav/OUT/Spectrum_Nphs_wmat1.75995007E-02_wcav5.86650008E-03_lamb3.00000003E-03.txt'
offset3 = 0
replot File3 u 1:2 w l lw 2 title '3 photon(s)'

#set key left top
