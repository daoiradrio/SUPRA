ir_for_exp_file = "final_ir_for_exp.dat"
ir_for_sim_file = "final_ir_for_sim.dat"
glob_min_ir_file = "final_glob_min_ir.dat"
exp_ir_file = "final_exp_ir.dat"
sim_ir_file = "final_sim_ir.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "dehydro_ir_spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set multiplot layout 3,1

set label 1 "{/Symbol s} = 0.9593" font "arial,10" at 1364,0.25
set label 2 "Scaling factor = 0.985" font "arial,10" at 1300,0.1
set key at 1220,0.3
set xlabel "wavenumber [1/cm]"
set ylabel "intensity" 
set xrange [1000:1700]
set yrange [-0.1:1.1]
plot \
ir_for_exp_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
exp_ir_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.9021" font "arial, 10" at 1357,0.25
set label 2 "Scaling factor = 0.97" font "arial,10" at 1300,0.1
set key at 1220,0.3
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1700]
set yrange [-0.1:1.1]
plot \
glob_min_ir_file u 1:2 w l lt 7 lw 2 title "computed global minimum", \
exp_ir_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.9696" font "arial, 10" at 1357,0.25
set label 2 "Scaling factor = 0.97" font "arial,10" at 1300,0.1
set key at 1220,0.3
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1700]
set yrange [-0.1:1.1]
plot \
ir_for_sim_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
sim_ir_file u 1:2 w l lt 8 lw 2 title "computed reference"

unset multiplot
