#!/usr/bin/gnuplot -persist
set terminal jpeg font arial 12 size 800,600
set output "test_g/DoS-L=16/graphs/15.jpg"
set grid x y
set xlabel "i"
set ylabel "G(i)"
plot "test_g/DoS-L=16/15.dat" using 1:3 title "landau-wang-16-iteration-15" with lines lt rgb "red"