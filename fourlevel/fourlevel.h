#ifndef FOURLEVEL_H
#define FOURLEVEL_H

#include <complex>
#include "/home/asudler/git/capstone/linalg/matrix.h"
#include "/home/asudler/git/capstone/spline/spline.h"

struct fourlevel_state
{
    // idk if this is allowed
    matrix<std::complex<double>> H2(double t);

    matrix<std::complex<double>> H, rho;
    double hbar, cap_gamma, ti, tf, cap_omega_plus, cap_omega_pi, 
           cap_omega_minus, cap_delta_B, cap_delta_pi, cap_delta_plus,
           cap_delta_upper; // should be const?

    fourlevel_state(); // default constructor
    fourlevel_state
    (
        matrix<std::complex<double>> H_init, 
        matrix<std::complex<double>> rho_init, 
        double hbar_init,
        double cap_gamma_init,
        double ti_init,
        double tf_init
    ); // parameterized constructor
    fourlevel_state(std::string inputfile); // ctr via ifstream
    ~fourlevel_state() {} // destructor

    /* member functions */
    std::pair<std::vector<double>, 
        std::vector<std::vector<std::complex<double>>>> solve(); 

// std::vector<cubic_spline<std::complex<double>>> solve();
// maybe come up with a better interface in the future
// so that we can easily do either the spline or the nonsplined rkf45

    private:
        matrix<std::complex<double>> H_default_();
        matrix<std::complex<double>> rho_default_();
};

//std::vector<std::complex<double>> master_eq(double, fourlevel_state);

#endif // FOURLEVEL_H
