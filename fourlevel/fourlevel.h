#ifndef FOURLEVEL_H
#define FOURLEVEL_H

#include <complex>
#include "/home/asudler/git/capstone/linalg/matrix.h"

struct fourlevel_state
{
    matrix<std::complex<double>> H, rho;
    double hbar, cap_gamma; // should be const?

    fourlevel_state(); // default constructor
    fourlevel_state
    (
        matrix<std::complex<double>> H_init, 
        matrix<std::complex<double>> rho_init, 
        double hbar_init,
        double cap_gamma_init
    ); // parameterized constructor
    fourlevel_state(std::string inputfile); // ctr via ifstream
    ~fourlevel_state() {} // destructor

    private:
        matrix<std::complex<double>> H_default_();
        matrix<std::complex<double>> rho_default_();
};

#endif // FOURLEVEL_H
