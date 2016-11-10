#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=6/graphs/17.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=6/17.dat" using 1:3 title "landau-wang-6-iteration-17" with lines lt rgb "red"