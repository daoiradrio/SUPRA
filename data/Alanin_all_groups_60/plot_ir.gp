datafile = "weighted_ir.dat"

set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 800, 600

set output "ir.png"
#set terminal pdf color font "Times-Roman, 10"
#set encoding iso_8859_1
#set key font "Times-Roman, 9"

set xlabel "wavenumber [1/cm]"
set ylabel "intensity"

plot datafile u 1:2 w l lt 7 lw 2 notitle
