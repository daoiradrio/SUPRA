ir_file = "normalized_scaled_glob_min_ir.dat"
ref_ir_file = "normalized_sim_ir.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "testplot.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set label 1 "{/Symbol s} = 0.2400" at 1600,1.3 font "arial, 13"

set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [950:1750]
set yrange [-0.1:1.1]
plot \
ir_file u 1:2 w l lt 7 lw 2 notitle, \
ref_ir_file u 1:2 w l lt 8 lw 2 notitle
