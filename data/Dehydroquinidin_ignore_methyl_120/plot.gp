#ir_file = "normalized_inverted_sliced_ir.dat"
#exp_ir_file = "normalized_ir_exp.dat"
#sim_ir_file = "sim_ir.dat"
vcd_file = "normalized_sliced_vcd.dat"
#exp_vcd_file = "exp_vcd.dat"
sim_vcd_file = "sim_vcd.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [950:1750]
set yrange [-1.5:1.5]
plot \
vcd_file u 1:2 w l lt 7 lw 2 notitle, \
sim_vcd_file u 1:2 w l lt 3 lw 2 notitle
#exp_vcd_file u 1:2 w l lt 2 lw 2 notitle, \
