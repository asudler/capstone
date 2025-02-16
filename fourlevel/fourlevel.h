#ifndef FOURLEVEL_H
#define FOURLEVEL_H

#include <complex>
#include <functional>
#include "/home/asudler/git/capstone/linalg/matrix/matrix.h"
#include "/home/asudler/git/capstone/spline/spline.h"

struct filenames
{
    std::string beams, rho, rho_log, spatial_log, spatial_rho_xii_log,
        spatial_rho_xif_log, spatial_omega_base, polaritons_base;
    int omega_print;

    filenames();
    filenames(std::string inputfile);
    ~filenames() {}
};

/* boundary_conditions:
 * it contains information about the atomic medium */
struct boundary_conditions
{
    double xi_min, xi_max, mu_alpha, epsilon, tolerance;
    int nxi, nn;

    boundary_conditions(); // default constructor
    boundary_conditions(std::string inputfile); // ctr via ifstream
    ~boundary_conditions() {} // destructor
}; // boundary_conditions


struct fourlevel_state
{
    matrix<std::complex<double>> H(double), rotate(double), rho0;
    cubic_spline<std::complex<double>> cap_omega_plus_t, cap_omega_pi_t,
        cap_omega_minus_t;
    std::function<double(double)> theta, phi;
    std::vector<double> times;
    std::vector<std::vector<std::complex<double>>> solutions;
//    std::vector<cubic_spline<std::complex<double>>> solutions_spline;

    double hbar, cap_gamma, cap_omega_plus, cap_omega_pi, cap_omega_minus, 
           cap_delta_B, cap_delta_pi, cap_delta_plus, cap_delta_upper, ti,
           tf, dt, t_on_pi, t_off_pi, tau_pi, t_on1_pm, t_off1_pm, 
           t_on2_pm, t_off2_pm, tau_pm, t_B_on, g, chi_m, chi_p;
    int nt, const_dt, nn;

    fourlevel_state(); // default constructor
    fourlevel_state
    (
        matrix<std::complex<double>> rho0_init, 
        double hbar_init,
        double cap_gamma_init,
        double cap_omega_plus_init,
        double cap_omega_pi_init, 
        double cap_omega_minus_init,
        double cap_delta_B_init,
        double cap_delta_pi_init,
        double cap_delta_plus_init,
        double cap_delta_upper_init,
        cubic_spline<std::complex<double>> cap_omega_plus_t_init,
        cubic_spline<std::complex<double>> cap_omega_pi_t_init,
        cubic_spline<std::complex<double>> cap_omega_minus_t_init
    ); // parameterized constructor
    fourlevel_state(std::string inputfile); // ctr via ifstream
    ~fourlevel_state() {} // destructor

    /* member functions */
    //std::pair<std::vector<double>, 
    //    std::vector<std::vector<std::complex<double>>>> 
    //    solve();
    void solve(void);
    void print_beams(std::string);
    void print_rabi_couplings(std::string);
    void print_rho(std::string);
    void print_rho_log(std::string);
    void print_polaritons(std::string);
    void print_polaritons_log(std::string);

    private:
        matrix<std::complex<double>> rho0_default_();
        cubic_spline<std::complex<double>> cap_omega_spline_default_();
        std::vector<std::complex<double>> master_eq_
        (
            double, 
            std::vector<std::complex<double>>
        );
}; // fourlevel_state

#endif // FOURLEVEL_H
