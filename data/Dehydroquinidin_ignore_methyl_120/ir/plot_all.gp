ir_file = "scaled_aligned_ir.dat"
glob_min_ir_file = "scaled_aligned_glob_min_ir.dat"
exp_ir_file = "aligned_exp_ir.dat"
sim_ir_file = "aligned_sim_ir.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "dehydro_ir_spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set multiplot layout 3,1

set label 1 "{/Symbol s} = 0.8906" font "arial,10" at 1320,0.05
set key at 1220,0.3
set xlabel "wavenumber [1/cm]"
set ylabel "intensity" 
set xrange [1000:1700]
set yrange [-0.1:1.1]
plot \
ir_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
exp_ir_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.8185" font "arial, 10" at 1320,0.05
set key at 1220,0.3
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1700]
set yrange [-0.1:1.1]
plot \
glob_min_ir_file u 1:2 w l lt 7 lw 2 title "computed global minimum", \
exp_ir_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.9016" font "arial, 10" at 1320,0.05
set key at 1220,0.3
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1700]
set yrange [-0.1:1.1]
plot \
ir_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
sim_ir_file u 1:2 w l lt 8 lw 2 title "computed reference"

unset multiplot
