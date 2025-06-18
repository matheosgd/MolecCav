reset
set term qt font "Times, 12"
unset grid
show grid

unset title
set xlabel 'Transition energy [cm-1]'
set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]

unset xrange
unset yrange
#min = 5.8304771453041417E-003 - 1E-4
#max = 5.9023030045979455E-003 + 1E-4
#Conv = 219474.6
#set xrange [min*Conv:max*Conv]
#set yrange [0.0249:0.049]

#set key left top

File1 = '/home/segaud/MolecCav/OUT/Spectrum_Nphs_wmat1.88565999E-02_wcav1.88565999E-02_lamb3.00000003E-03.txt'
offset1 = 0
  plot File1 u 1:2 w l lw 2 title '1 photon(s)'

File2 = '/home/segaud/MolecCav/OUT/Spectrum_Nphs_wmat1.88565999E-02_wcav3.77131999E-02_lamb3.00000003E-03.txt'
offset2 = 0
replot File2 u 1:2 w l lw 5 title '2 photon(s)'

File3 = '/home/segaud/MolecCav/OUT/Spectrum_Nphs_wmat1.88565999E-02_wcav5.65697998E-02_lamb3.00000003E-03.txt'
offset3 = 0
replot File3 u 1:2 w l lw 2 title '3 photon(s)'

#set key left top
