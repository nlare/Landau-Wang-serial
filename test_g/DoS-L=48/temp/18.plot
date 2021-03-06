#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=48/graphs/18.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=48/18.dat" using 1:3 title "landau-wang-48-iteration-18" with lines lt rgb "red"