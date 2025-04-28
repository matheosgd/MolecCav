reset
set term qt font "Times, 12"
set grid
show grid

set title  'Intensities = f(coupling\_strength) for resonant case (w = w_{HF}) [a.u.]'
set xlabel 'Transition energy [Ha]'
set ylabel 'Arbitrary units' #Coupling strength [a.u.]

unset xrange
unset yrange
min = 5.7455455693706022E-003 - 1E-4
max = 5.9850104937492297E-003 + 1E-4
set xrange [min:max]
#set yrange [0.0249:0.049]
#set key left top

Conv = 21947.46 # 1Ha = Conv.cm-1
Gam = 0.1/Conv  # 10cm-1
set samples 800
show samples
L(x, x_0, Gam) = ( Gam/(2*pi) ) / ( ((Gam**2)/4) + (x-x_0)**2 )

offset = 0 #0.000*1E12
  plot ( 0.0488532887*L(x, 5.8665000000000002E-003, Gam) + 0.0000000000*L(x, 5.8665000000000002E-003, Gam) + offset) w l lw 2 t '\lambda = 0'

offset = 0 #0.002*1E12
replot ( 0.0245263300*L(x, 5.8426073486058791E-003, Gam) + 0.0243269587*L(x, 5.8904903575647703E-003, Gam) + offset) w l lw 2 t '\lambda = 0.002'

offset = 0 #0.004*1E12
replot ( 0.0246260106*L(x, 5.8188124009415547E-003, Gam) + 0.0242272781*L(x, 5.9145784188593372E-003, Gam) + offset) w l lw 2 t '\lambda = 0.004'

offset = 0 #0.006*1E12
replot ( 0.0247256813*L(x, 5.7951151496850734E-003, Gam) + 0.0241276074*L(x, 5.9387641765617505E-003, Gam) + offset) w l lw 2 t '\lambda = 0.006'

offset = 0 #0.008*1E12
replot ( 0.0248253371*L(x, 5.7715155826352210E-003, Gam) + 0.0240279516*L(x, 5.9630476184707893E-003, Gam) + offset) w l lw 2 t '\lambda = 0.008'

offset = 0 #0.010*1E12
replot ( 0.0249249729*L(x, 5.7480136827145539E-003, Gam) + 0.0239283158*L(x, 5.9874287275090135E-003, Gam) + offset) w l lw 2 t '\lambda = 0.010'

#set key left top
