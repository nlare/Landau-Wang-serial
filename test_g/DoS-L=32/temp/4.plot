#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=32/graphs/4.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=32/4.dat" using 1:3 title "landau-wang-32-iteration-4" with lines lt rgb "red"