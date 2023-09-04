vcd_for_exp_file = "final_vcd_for_exp.dat"
vcd_for_sim_file = "final_vcd_for_sim.dat"
glob_min_vcd_file = "final_glob_min_vcd.dat"
exp_vcd_file = "final_exp_vcd.dat"
sim_vcd_file = "final_sim_vcd.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "dehydro_vcd_spectra.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set multiplot layout 3,1

set label 1 "{/Symbol s} = 0.2513" font "arial,10" at 1622,-0.65
set label 2 "Scaling factor = 0.985" font "arial,10" at 1558,-1
set key on
set xlabel "wavenumber [1/cm]"
set ylabel "intensity" 
set xrange [1000:1700]
set yrange [-1.5:1.5]
plot \
vcd_for_exp_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
exp_vcd_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.0437" font "arial, 10" at 1622,-0.65
set label 2 "Scaling factor = 0.98" font "arial,10" at 1565,-1
set key on
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1700]
set yrange [-1.5:1.5]
plot \
glob_min_vcd_file u 1:2 w l lt 7 lw 2 title "computed global minimum", \
exp_vcd_file u 1:2 w l lt 8 lw 2 title "experimental reference"

set label 1 "{/Symbol s} = 0.4178" font "arial, 10" at 1622,-0.65
set label 2 "Scaling factor = 0.98" font "arial,10" at 1565,-1
set key inside
set xlabel "wavenumber [1/cm]"
set ylabel "intensity"
set xrange [1000:1700]
set yrange [-1.5:1.5]
plot \
vcd_for_sim_file u 1:2 w l lt 7 lw 2 title "computed ensemble", \
sim_vcd_file u 1:2 w l lt 8 lw 2 title "computed reference"

unset multiplot
