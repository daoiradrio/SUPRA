ir_file = "normalized_inverted_sliced_ir.dat"
vcd_file = "normalized_sliced_vcd.dat"
exp_ir_file = "normalized_ir_exp.dat"
exp_vcd_file = "normalized_vcd_exp.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set multiplot layout 2,1

set label 1 'IR' at graph 0.92,0.9 font ',8'
set xlabel "wavenumber [1/cm]"
set ylabel "intensity" 
plot \
ir_file u 1:2 w l lt 7 lw 2 notitle, \
exp_ir_file u 1:2 w l lt 2 lw 2 notitle

set label 1 'VCD' at graph 0.92,0.9 font ',8'
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [200:2000]
set yrange [-2.5:1.5]
plot \
vcd_file u 1:2 w l lt 7 lw 2 notitle, \
exp_vcd_file u 1:2 w l lt 2 lw 2 notitle

unset multiplot
