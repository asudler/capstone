#ifndef FOURLEVEL_H
#define FOURLEVEL_H

#include <complex>
#include "/home/asudler/git/capstone/linalg/matrix.h"
#include "/home/asudler/git/capstone/spline/spline.h"

struct fourlevel_state
{
    matrix<std::complex<double>> H(double t), rho;
    double hbar, cap_gamma, cap_omega_plus, cap_omega_pi, cap_omega_minus, 
           cap_delta_B, cap_delta_pi, cap_delta_plus, cap_delta_upper, ti, tf,
           t_on, t_off, tau; 
    // should be const?

    fourlevel_state(); // default constructor
    fourlevel_state
    (
        matrix<std::complex<double>> rho_init, 
        double hbar_init,
        double cap_gamma_init,
        double cap_omega_plus_init,
        double cap_omega_pi_init, 
        double cap_omega_minus_init,
        double cap_delta_B_init,
        double cap_delta_pi_init,
        double cap_delta_plus_init,
        double cap_delta_upper_init,
        double ti_init,
        double tf_init,
        double t_on_init,
        double t_off_init,
        double tau_init
    ); // parameterized constructor
    fourlevel_state(std::string inputfile); // ctr via ifstream
    ~fourlevel_state() {} // destructor

    /* member functions */
    std::pair<std::vector<double>, 
        std::vector<std::vector<std::complex<double>>>> solve(); 
    // future development should include splined solution

    private:
        matrix<std::complex<double>> rho_default_();
};

#endif // FOURLEVEL_H
