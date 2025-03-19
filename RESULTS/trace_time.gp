reset
set term qt font "Times, 15"
#set term png font "Times, 15"
#set output "/home/msegaud/MolecCav/RESULTS/time.png"
set grid
show grid

#set datafile separator ";"
set xlabel 'Nb_C = Nb_M'
set ylabel 'Time (cpu)'
#set zlabel 'z'
#set title  'Ratio (Absorbance du SQ pour diff. [H_{2}O_{2}]) / (Abs(SQ seul))'
unset xrange
unset yrange
#unset zrange
#set yrange [1.599960484:1.599960489]
#set xrange [0:55]
#set zrange [-10:10]

#offset1 = 0
#offset2 = 0
#offset3 = 0
#scale   = (1/0.225)

#EX :
#DISPLAY THE TEST PLOT (in the gnuplot terminal)
#   gnuplot> test
#PLOT CLASSIQUE + SUPERPOSER GRAPHS + CHOISIR COULEUR
#     plot '/home/msegaud/outils/CSV/PB3_1mM_200-900_06.csv' skip 19 using 1:($2*scale+offset1) with lines lw 2 lc 'blue' title '[H_{2}O_{2}] = 1mM'
#   replot '/home/msegaud/outils/CSV/PB3_1mM_200-900_06.csv' skip 19 u 1:3 w l lw 2 t '[H_{2}O_{2}] = 1mM'
#Rq: /!\ La numérotation des colonnes COMMENCE à 1 !!!!!!!!!! /!\
#Rq: if only one column is provided after the using, these values will be plotted in ordinates and the abscissa will be the number of the point (=the line most of the time).
#PLOT CLASSIQUE D'UNE FONCTION DIRECTEMENT DEFINIE
#   plot cos(x)
#Rq : log : base e; log10 : base 10
#PLOT CLASSIQUE D'UNE FONCTION DEFINIE AVANT
#   deltaH            = 1000
#   deltaS            = 5
#   gamma_interaction = 250
#   R                 = 8.314
#   f(x)              = (deltaH + gamma_interaction*(1-2*x)) / (R*log((1-x)/x)+deltaS)
#   plot f(x) 
#SUPERPOSE PLS COURBES EX PLUSIEURS FICHIERS
#   offset1 = 0
#   offset2 = 0
#   offset3 = 0
#   scale   = (1/0.108)
#     plot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($2*scale+offset1) w l lw 2 t 'x' 
#   replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($3*scale+offset2) w l lw 2 t 'x'
#   replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($4*scale+offset3) w l lw 2 t 'x'
#TRACER MAIS SANS AFFICHER POINTS/LIGNE
#   plot ... lt rgb 'white' (pour lignes) (peut être que ça marche aussi ac lc 'white' en fait)
#   plot ... pointtype -1   (pour points)
#PLOT POINTS
#   plot '/home/msegaud/outils/CSV/VoltaCoucheFineH2O2.csv' u 1:2 w points pointtype 7 pointsize 1.5 lc 'blue' t '[H_{2}O_{2}] = 0'
#PLOT LINES AND POINTS
#   plot '/home/msegaud/outils/CSV/VoltaCoucheFineH2O2.csv' w linespoints lc 'cyan' lw 2 pointtype 5
#PLOT ERRORBARS ET SANS TITRE
#   plot '/home/msegaud/outils/CSV/ratio_couche_fine.csv' u 1:2:3 with <x or y>errorbars pointtype 7 pointsize 1 notitle #t '[H_{2}O_{2}] = 0.5mM'
#Rq: to get lines plotted between the data points, plot the data twice, once with errorbars and once with lines
#Rq: les arrorbars sont parametrées comme si elles étaient des lignes : modifier lw, lc, dt (dashtype)
#Rq: peut ne tracer que l'errorbar (sans les points centraux) en enlevant simplement pointtype et pointsize de la ligne
#Rq: plutôt que de parametrer les caracteristiques de la ligne dans le plot peut le faire par set errorbars : 
#   set errorbars <small|large> <line properties>
#Rq: small enlève les petits traits aux extrémités
#PLOT AVEC VARIABLE DANS LA LEGENDE
#   plot '/home/msegaud/outils/CSV/A_PB_206nm.csv' u 1:($2*scale) w points pointtype 7 pointsize 1.5 notitle, f(x) t sprintf("R^2 = %g", STATS_correlation) #notitle
#SUPERPOSER PLOTS 3D EX PLUSIEURS FICHIERS
#    splot   '/home/msegaud/outils/CSV/PI_complex_catal_opt' u 1:2:3 w points pointtype 7 pointsize 1.5 lc 'blue'  t 'catal'
#   replot  '/home/msegaud/outils/CSV/PI_complex_eth_opt'   u 1:2:3 w points pointtype 7 pointsize 1.5 lc 'green' t 'ethylene'
#CHOISIR LES POINTS A TRACER
#                                                      every <pt_inc>:<b_inc>:<st_pt>:<st_bl>:<end_pt>:<end_bl>
#    splot   '/home/msegaud/outils/CSV/PI_complex_opt' every         :       :0      :       :0                 u 1:2:3 w points pointtype 7 pointsize 1.5 lc 'blue'   t 'Ti'
#   replot   '/home/msegaud/outils/CSV/PI_complex_opt' every         :       :1      :       :1                 u 1:2:3 w points pointtype 7 pointsize 1.5 lc 'black'  t 'C (Me)'
#Rq: /!\ La numérotation des points/blocs COMMENCE à 0 !!!!!!!!!! /!\
#PM3D
#   set pm3d at bs explicit #s:surface; b:projection at the bottom; t:projection at the top
#   splot '/home/msegaud/Internship_CBBC/olefin_polymn/orbitals/PI_complex_opt_6-31Gdp_GenWIa.mo10a.txt' every ::13  with pm3d                                         notitle
#CHANGER SYST COO EN SPLOT
#   set mapping <cartesian (default) | spherical | cylindrical>
#Rq: spherical : col 1 : theta; col 2 : phi; col 3 : radius (default = 1)
#Rq: cylindrical : col 1 : theta; col 2 : z; col 3 : radius (default = 1)
#Rq: penser à changer le nom des var pour rendre naturel l'écriture si ecrit des fonctions, e.g. pour cylindrical :
#   set dummy theta, z, r
#METTRE AUTRE ABSCISSE
#   set xtics
#     plot '/home/msegaud/outils/CSV/reaction_path' every :::0::0 u 1:(($2-offset1)*scale):xticlabels(3) w points pointtype 5 pointsize 1.5 t 'BP86'
#   replot '/home/msegaud/outils/CSV/reaction_path' every :::1::1 u 1:(($2-offset2)*scale):xticlabels(3) w points pointtype 9 pointsize 1.5 t 'RMP2'
#Rq: (le every ici est juste par ce qu'on avait deux blocs mais aucun lien ac xtics)
#Rq: data file = col 1 : abscisse; col 2 : valeur; col 3 : nouveaux noms des points de l'absc. EX : 
#1 -968.9437494066269 R
#2 -968.99563898384 Pi
#...
#6 -900 
#REGRESSION LINEAIRE/FIT AVEC FONCTION
#   stats '/home/msegaud/outils/CSV/A_PB_206nm.csv' using 1:2 #Cree les valeurs stats (STATS_correlation)
#   f(x) = a*x + b
#   set fit logfile "ratio_PB_fit" #Indique dans quel fichier ecrire les resultats du fit
#   fit f(x) '/home/msegaud/outils/CSV/A_PB_206nm.csv' u 1:($2*scale) via a,b #le tracé
#AJOUTER UNE DROITE VERTICALE DANS LE PLOT
#   set arrow from 206, graph 0 to 206, graph 1 nohead #trace droite verticale à abscisse 206 nm
#AFFICHER OU NON LA GRILLE
#   set grid
#   show grid
#Rq: set grid <peramètres> supporte plusieurs arguments optionels et permet en fait de la parametrer.
#DEPLACER LA LEGENDE
#   set key right bottom #deplace legende (=key) en bas à droite /!\ penser à remettre en haut à la fin (cf derniere ligne)
#DEFINIR UN STYLE DE PLOT
#   set style line 1 lt rgb 'cyan' lw 3 dt 2 pt 5
#   plot ... with linespoints ls 1
#Rq: 1 le nom associé au style par lequel il sera invoqué; lt linetype (?); lw linewidth; dt dashtype; pt pointtype
#Rq: ls linestyle
#REDIRIGER LA SORTIE DANS UN FICHIER (i.e. sauver la figure sans l'afficher)
#   set term png font "Times, 15"
#   set output "/home/msegaud/MolecCav/RESULTS/convergence_ZPE.png"
#/!\/!\/!\ les replot ne marchent pas avec une redirection /!\/!\/!\
#Rq: Ceci pour une figure .png, si autre (gif par ex.), se réferer au terminal correspondant plutôt que png. 
#MEILLEUR TERMINAL POSSIBLE
#   set term dumb size <nb of xchar>, <nb of ychars> (optional) aspect <ratio hlenght/vlength> (optional) <mono|ansirgb> (optional)
#Rq: mono : noir et blanc (default); ansirgb sequence de couleurs classiques de gnuplot. Il existe d'autres options (ansi et ansi256).
#Rq: /!\ les points sur ligne abscisses ne seront pas tracés : mieux vaut décaler l'origine du repère.
#Rq: si deux points sont trop il va les paniquer et les mélanger : juste augmenter la taille du graph avec size ! (default: 79,24)
#TRACER EN ECHELLE LOG 
#Rq: deux options : 1) tracer u 1:(log($2)) 2) changer parametres d'échelle
#   plot '/home/msegaud/MolecCav/RESULTS/time_test_02.txt' every :::0::0 using 1:(log10($2)) with points pointtype 7 pointsize 1.5 t 'time\_dense'
#ou
#   set logscale y
#   plot '/home/msegaud/MolecCav/RESULTS/time_test_02.txt' every :::0::0 using 1:2 with points pointtype 7 pointsize 1.5 t 'time\_dense'


#set term dumb size 100, 30 feed ansirgb
set xrange [0:55]
#set term png font "Times, 15"
#set output "/home/msegaud/MolecCav/RESULTS/time.png"
#set logscale y
#set logscale x
set key left top #deplace legende (=key) en bas à droite /!\ penser à remettre en haut à la fin (cf derniere ligne)
  plot '/home/msegaud/MolecCav/RESULTS/time_test_02.txt' every :::0::0 using 1:2 with points pointtype 7 pointsize 1.5 t 'time\_opt'
replot '/home/msegaud/MolecCav/RESULTS/time_test_02.txt' every :::1::1 using 1:2 with points pointtype 7 pointsize 1.5 t 'time\_dense'
set key right top #reinitialise position legende après les tracé


#set key right top #reinitialise position legende après les tracé
