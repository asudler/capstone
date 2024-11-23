unset datafile separator
p 'out.rho.dat' u 1:2 title "rho11", '' u 1:12 title "rho22", \
'' u 1:22 title "rho33", '' u 1:32 title "rho44", \
'../../../../doerte/basic-four-level/densitymatrix_11.dat' u 1:3 w l title \
"doerte rho11", \
'../../../../doerte/basic-four-level/densitymatrix_22.dat' u 1:3 w l title \
"doerte rho22", \
'../../../../doerte/basic-four-level/densitymatrix_33.dat' u 1:3 w l title \
"doerte rho33", \
'../../../../doerte/basic-four-level/densitymatrix_77.dat' u 1:3 w l title \
"doerte rho77"
