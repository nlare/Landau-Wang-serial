#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=64/graphs/6.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=64/6.dat" using 1:3 title "landau-wang-64-iteration-6" with lines lt rgb "red"