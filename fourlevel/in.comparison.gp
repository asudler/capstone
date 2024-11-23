set grid
set xlabel "time"
set ylabel "|\Omega_{pi, nxi}|^2"

p [2.5:4.5][0:0.00025] 'test-grid/out.plot_omega_pi_nxi.dat' u 1:(($2)**2 + ($3)**2) w l title "non-const dt", 'test-grid-const-dt/out.plot_omega_pi_nxi.dat' u 1:(($2)**2 + ($3)**2) w l title "const dt", '< paste ../../../doerte/in-medium/fort.30010 ../../../doerte/in-medium/fort.30011' u 1:(($2)**2 + ($4)**2) w l title "doerte" 
