reset
set term qt font "Times, 12"
set grid
show grid

set title  'Intensities = f(detuning) for small coupling (\lambda = 4e-4 a.u.; w_M = w_{HF}) [a.u.]
set xlabel 'Transition energy [Ha]'
set ylabel 'Arbitrary units [a.u.]'

unset xrange
unset yrange
min = 5.4610394673397069E-003 - 1E-5
max = 6.2723303490619256E-003 + 1E-5
set xrange [min:max]
#set yrange [0.0249:0.049]
#set key left top

Conv = 21947.46 # 1Ha = Conv.cm-1
Gam = 0.3/Conv  # 30cm-1
set samples 999
show samples
L(x, x_0, Gam) = ( Gam/(2*pi) ) / ( ((Gam**2)/4) + (x-x_0)**2 )

offset = 0 #0*1E12
  plot ( 0.0246284919*L(x, 5.8184199664414512E-003, Gam) + 0.0242296797*L(x, 5.9141891746978970E-003, Gam) + offset) w l lw 2 t 'Detuning = 0'

offset = 0 #.0002*1E12
replot ( 0.0464792636*L(x, 5.8550758224022222E-003, Gam) + 0.0023789621*L(x, 6.0775267678491625E-003, Gam) + offset) w l lw 2 t 'Detuning = .0002'

offset = 0 #-.0002*1E12
replot ( 0.0024110328*L(x, 5.6557875323548475E-003, Gam) + 0.0464470819*L(x, 5.8768283868815456E-003, Gam) + offset) w l lw 2 t 'Detuning = -.0002'
#
set key left top
