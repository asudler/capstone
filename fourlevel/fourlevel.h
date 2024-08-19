#ifndef FOURLEVEL_H
#define FOURLEVEL_H

#include <complex>
#include "/home/asudler/git/capstone/linalg/matrix.h"

struct fourlevel_state
{
    matrix<std::complex<double>> H, rho;
    double hbar, cap_gamma, ti, tf; // should be const?

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

    private:
        matrix<std::complex<double>> H_default_();
        matrix<std::complex<double>> rho_default_();
};

//std::vector<std::complex<double>> master_eq(double, fourlevel_state);

#endif // FOURLEVEL_H
