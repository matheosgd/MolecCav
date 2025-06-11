reset
set term qt font "Times, 12"                                                           
set grid
show grid

File1 = '/home/segaud/MolecCav/OUT/Spec_1p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
File2 = '/home/segaud/MolecCav/OUT/Spec_2p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'
File3 = '/home/segaud/MolecCav/OUT/Spec_3p1D_wmat5.86650008E-03_wcav5.86650008E-03_lamb3.00000003E-03.txt'

Conv = 219474.6 # 1Ha = Conv.cm-1
#Gam = 1.0/Conv # 10cm-1 (line useful if plot in Ha)
Gam = 1.0       # 10cm-1 (line useful if plot in cm-1)
L(x, x_0, Gam) = ( Gam/(2*pi) ) / ( ((Gam**2)/4) + (x-x_0*Conv)**2 )
offset = 0 #coupling_strength*1E12

unset title #set title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'
set xlabel 'Transition energy [Ha]'
set ylabel 'Transition intensity [Arbitrary units]' #Coupling strength [a.u.]

unset xrange
unset yrange
min = 5.8154906356422710E-003 - 1E-4
max = 5.9170696435614355E-003 + 1E-4
set xrange [min:max]
#set xrange [min*Conv:max*Conv]
#set yrange [0.0249:0.049]

#set key left top
set samples 900
show samples

  plot File1 u 1:2 w point lw 2 title '1p1D'
replot File2 u 1:2 w point lw 2 title '2p1D'
replot File3 u 1:2 w point lw 2 title '3p1D'
#  plot ( ${GSto1}*L(x, ${Enrgy1}, Gam) + ${GSto2}*L(x, ${Enrgy2}, Gam) + offset) w l lw 2 t '\lambda = 0'
#replot ( ${GSto1}*L(x, ${Enrgy1}, Gam) + ${GSto2}*L(x, ${Enrgy2}, Gam) + offset) w l lw 2 t '\lambda = $coupling_strength'

#set key left top