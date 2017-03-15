#!/bin/bash
gnuplot test_g/DoS-L=32/temp/*.plot
convert -delay 100 -loop 0 test_g/DoS-L=32/graphs/{1..20}.jpg animate-DoS-L=32.gif
