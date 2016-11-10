#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=6/graphs/20.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=6/20.dat" using 1:3 title "landau-wang-6-iteration-20" with lines lt rgb "red"