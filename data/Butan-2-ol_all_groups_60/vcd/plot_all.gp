vcd_file = "aligned_scaled_vcd.dat"
glob_min_vcd_file = "aligned_scaled_glob_min_vcd.dat"
exp_vcd_file = "aligned_exp_vcd.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "butanol_vcd_spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set multiplot layout 2,1

set label 1 "{/Symbol s} = 0.4931" font "arial,10" at 1312,1
set label 2 "Scaling factor = 0.985" font "arial,10" at 1277,0.8
set key on
set xlabel "wavenumber [1/cm]"
set ylabel "intensity" 
set xrange [1000:1375]
set yrange [-1.1:1.7]
plot \
vcd_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
exp_vcd_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.2669" font "arial, 10" at 1312,1
set label 2 "Scaling factor = 0.98" font "arial,10" at 1281,0.8
set key on
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1375]
set yrange [-1.1:1.7]
plot \
glob_min_vcd_file u 1:2 w l lt 7 lw 2 title "computed global minimum", \
exp_vcd_file u 1:2 w l lt 8 lw 2 title "experimental reference"

unset multiplot
