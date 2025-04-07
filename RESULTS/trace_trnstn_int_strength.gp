reset
set term qt font "Times, 12"
set grid
show grid

set title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'
set xlabel 'Transition energy [Ha]'
set ylabel 'Coupling strength [a.u.]'

unset xrange
unset yrange
#set xrange [4.4E-003:7.0E-003]
#set yrange [0.0249:0.049]
#set key left top

Conv = 21947.46 # 1Ha = Conv.cm-1
Gam = 0.1/Conv  # 10cm-1
L(x, x_0, Gam) = ( Gam/(2*pi) ) / ( ((Gam**2)/4) + (x-x_0)**2 )
offset = 0#0.00*1E12

  plot ( 0.0488532887*L(x, 5.8665000000000002E-003, Gam) + 0.0000000000*L(x, 5.8665000000000002E-003, Gam) + offset) w l lw 2 t '\lambda = 0'

min = 5.8665000000000002E-003 - 1E-4
max = 5.8665000000000002E-003 + 1E-4
set xrange [min:max]
offset = 0#0.002*1E12
replot ( 0.0245269452*L(x, 5.8425094418358900E-003, Gam) + 0.0243275640*L(x, 5.8903928495522562E-003, Gam) + offset) w l lw 2 t '\lambda = 0.002'

min = 5.8425094418358900E-003 - 1E-4
max = 5.8903928495522562E-003 + 1E-4
set xrange [min:max]
offset = 0#0.004*1E12
replot ( 0.0246284919*L(x, 5.8184199664414512E-003, Gam) + 0.0242296797*L(x, 5.9141891746978970E-003, Gam) + offset) w l lw 2 t '\lambda = 0.004'

min = 5.8184199664414512E-003 - 1E-4
max = 5.9141891746978970E-003 + 1E-4
set xrange [min:max]
offset = 0#0.006*1E12
replot ( 0.0247313104*L(x, 5.7942303400734278E-003, Gam) + 0.0241329674*L(x, 5.9378901359128993E-003, Gam) + offset) w l lw 2 t '\lambda = 0.006'

min = 5.7942303400734278E-003 - 1E-4
max = 5.9378901359128993E-003 + 1E-4
set xrange [min:max]
offset = 0#0.008*1E12
replot ( 0.0248354275*L(x, 5.7699393031270417E-003, Gam) + 0.0240374041*L(x, 5.9614968706047213E-003, Gam) + offset) w l lw 2 t '\lambda = 0.008'

min = 5.7699393031270417E-003 - 1E-4
max = 5.9614968706047213E-003 + 1E-4
set xrange [min:max]
offset = 0#0.01*1E12
replot ( 0.0249408707*L(x, 5.7455455693706022E-003, Gam) + 0.0239429671*L(x, 5.9850104937492297E-003, Gam) + offset) w l lw 2 t '\lambda = 0.01'

min = 5.7455455693706022E-003 - 1E-4
max = 5.9850104937492297E-003 + 1E-4
set xrange [min:max]
#
set key left top
