set terminal pngcairo size 1024, 768 background "white"
set out outfile
set xlabel "x"
set ylabel "y"
set key outside right
set datafile separator ","

df(x)=exp(-x/10)*(-0.1*cos(x)-sin(x))
f(x)=exp(-x/10)*cos(x)
F(x)=-(10./101.)*exp(-x/10)*(cos(x)-10*sin(x))+0.1
set title "Cubic Spline (imag part), f(x) = exp(-x/10)*cos(x)"
plot f(x) title "function", \
df(x) title "derivative (analytical)", \
F(x) title "antiderivative (analytical)", \
data skip 1 u 1:3 w p pointsize 2 title "raw data", \
interpolation skip 1 u 1:3 w p title "cubic spline", \
interpolation skip 1 u 1:5 w p title "spline differentiation", \
interpolation skip 1 u 1:7 w p title "spline integration"

