vcd_file = "final_vcd_for_sim.dat"
exp_vcd_file = "final_sim_vcd.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "dehydro_spectra_own_exp.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set label 1 "{/Symbol s} = 0.1114" at 1600,1.3 font "arial, 13"

set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [950:1750]
set yrange [-1.5:1.5]
plot \
vcd_file u 1:2 w l lt 7 lw 2 notitle, \
exp_vcd_file u 1:2 w l lt 8 lw 2 notitle
