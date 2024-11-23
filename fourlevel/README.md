IMPORTANT!!!!! WE HARDWIRED A CHANGE TO THE RKF45 TOLERANCE IN FOURLEVEL.CC!

2024-09-21. Some comments on convergence.

Increasing the damping factor used for the medium boundary problem finite differencing from 0.7 --> 0.9 makes the calculation run much, much slower. Around 31 Gb of data were generated from a test run we performed earlier, so the directory has been deleted.

2024-10-14. More comments on convergence.

The following summarizes how the err and tol parameters on the RKF45 subroutine affect the convergence of the spatial grid results (obtained via iterative finite differencing). Note that epsilon (i.e. the damping factor on how much of the old pi-beam profile grid should be used in determining the new pi-beam profile grid) is 0.7 here. In what follows, 'norm' means the magnitude of all elements of the density matrix summed over all space and time points.

We have deleted the directories containing these calculations to save space. Ideally the comments below should contain enough description to get an idea of the best error and tolerance settings.

(err, tol) = (1e-3, 1e-3)
* the difference in norms between iterations begins around ~1e3 and oscillates rapidly
** these oscillations decrease in magnitude from it. 0 to it. 500
** a semi-linear (more aptly, linear with noise) decline begins around it. 500,
at which point the difference in norms is ~300
** the difference in the norms becomes less than 1e-3 after >=5750 its.
* the pi-beam profile grid agrees with Doerte's (and thus my const. dt) results
* unsure on the total timing of this run

(err, tol) = (1e-4, 1e-4)
* the difference in norms between iterations begins around ~1e2 and oscillates rapidly
** these oscillations decrease in magnitude from it. 0 to it. 500
** a semi-linear (more aptly, linear with noise) decline begins around it. 500, at which point the difference in norms is ~30
** the difference in the norms becomes less than 1e-3 after >=6600 its.
* the pi-beam profile grid agrees with Doerte's (and thus my const. dt) results
* unsure on the total timing of this run

(err, tol) = (1e-5, 1e-5)
* the difference in norms between iterations begins around ~1e1 and oscillates rapidly
** these oscillations decrease in magnitude from it. 0 to it. 500
** a semi-linear (more aptly, linear with noise) decline begins around it. 500,
at which point the difference in norms is ~3
** the difference in the norms becomes less than 1e-3 after >=6120 its.
* the pi-beam profile grid agrees with Doerte's (and thus my const. dt) results
* unsure on the total timing of this run

(err, tol) = (1e-6, 1e-6)
* !! something odd happened here. re-run at a later point !!
* !! the difference in norm values never converged !!

(err, tol) = (1e-9, 1e-9)
* the difference in the norms between iterations begins around ~0.25, and a nice linear decline is immediately observed (minimal noise)
** the convergence profile almost perfectly matches the const. dt results (although there is a bit of noise present in some areas, it is incredibly small)
** the difference in the norms becomes less than 1e-3 after >= 1720 its.
* the pi-beam profile grid agrees with Doerte's (and thus my const. dt) results
* unsure on the total timing of this run
