#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "graph/DensityStat_L=96/Hi.jpg"
set grid x y
set xlabel "i"
set ylabel "H(i)"
plot "results/DensityStat_L=96.dat" using 1:4 title "landau-wang-96" with lines lt rgb "red"