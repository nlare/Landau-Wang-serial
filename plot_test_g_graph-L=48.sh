#!/bin/bash
gnuplot test_g/DoS-L=48/temp/*.plot
convert -delay 100 -loop 0 test_g/DoS-L=48/graphs/{1..20}.jpg animate-DoS-L=48.gif
