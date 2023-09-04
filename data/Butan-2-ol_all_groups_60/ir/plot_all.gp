ir_file = "aligned_scaled_ir.dat"
glob_min_ir_file = "aligned_scaled_glob_min_ir.dat"
exp_ir_file = "aligned_exp_ir.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "butanol_ir_spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set multiplot layout 2,1

set label 1 "{/Symbol s} = 0.8367" font "arial,10" at 1302,-0.15
set label 2 "Scaling factor 0.985" font "arial,10" at 1270,-0.3
set key bottom right
set xlabel "wavenumber [1/cm]"
set ylabel "intensity" 
set xrange [950:1375]
set yrange [-0.8:1.1]
plot \
ir_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
exp_ir_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.7282" font "arial, 10" at 1302,-0.15
set label 2 "Scaling factor 0.985" font "arial,10" at 1270,-0.3
set key bottom right
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [950:1375]
set yrange [-0.8:1.1]
plot \
glob_min_ir_file u 1:2 w l lt 7 lw 2 title "computed global minimum", \
exp_ir_file u 1:2 w l lt 8 lw 2 title "experimental reference"

unset multiplot
