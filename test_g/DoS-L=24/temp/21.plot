#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=24/graphs/21.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=24/21.dat" using 1:3 title "landau-wang-24-iteration-21" with lines lt rgb "red"