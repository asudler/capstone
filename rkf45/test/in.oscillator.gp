set terminal pngcairo size 1024, 768 background "white"
set out outfile
set xlabel "t"
set key outside right

set title "Oscillator with Friction"
plot data skip 1 u 1:2 w lp title "x", \
'' skip 1 u 1:3 w lp title "dx/dt"

