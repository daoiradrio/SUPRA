vcd_file = "normalized_sliced_scaled_aligned_vcd.dat"
sim_vcd_file = "aligned_sim_vcd.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "dehydro_spectra_own_sim.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set label 1 "{/Symbol s} = 0.2400" at 1600,1.3 font "arial, 13"

set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [950:1750]
set yrange [-1.5:1.5]
plot \
vcd_file u 1:2 w l lt 7 lw 2 notitle, \
sim_vcd_file u 1:2 w l lt 6 lw 2 notitle
