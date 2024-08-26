#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include "fourlevel.h"
#include "/home/asudler/git/capstone/linalg/matrix.h"
#include "/home/asudler/git/capstone/misctools/misctools.h"
#include "/home/asudler/git/capstone/rkf45/rkf45.h"

using namespace std::complex_literals;
//double t_on = 2.3, t_off = 3.5, tau = 0.1;
/*double cap_omega_plus = 34.5575191895, cap_omega_pi = 15.7079632679,
       cap_omega_minus = 34.5575191895, cap_delta_B = 6.28318530718,
       cap_delta_pi = 12.5663706144, cap_delta_plus = 12.5663706144,
       cap_delta_upper = 0.;*/

// (1) initialize a state from input file parameters
// (2) solve the dynamics of this state
// (3) profit

matrix<std::complex<double>> fourlevel_state::H_default_()
{
    matrix<std::complex<double>> result(4,4);
    result(3,3) = 1.; // excited state energy = 1
    return result;
} // default Hamiltonian

matrix<std::complex<double>> fourlevel_state::rho_default_()
{
    matrix<std::complex<double>> result(4,4);
    result(0,0) = 1.; // initial state = ground
    return result;
} // default density matrix

fourlevel_state::fourlevel_state() : // default constructor
    H(H_default_()), rho(rho_default_()), hbar(1.), cap_gamma(1.),
    ti(0.), tf(1.) {}

fourlevel_state::fourlevel_state
(
    matrix<std::complex<double>> H_init,  
    matrix<std::complex<double>> rho_init, 
    double hbar_init,
    double cap_gamma_init,
    double ti_init,
    double tf_init
) : // parameterized constructor
    H(H_init), rho(rho_init), hbar(hbar_init), cap_gamma(cap_gamma_init),
    ti(ti_init), tf(tf_init) {}

fourlevel_state::fourlevel_state(std::string inputfile) : fourlevel_state()
{
    std::string rho_re_file, rho_im_file;
    
    std::ifstream istrm(inputfile, std::ios::binary);
    if(!istrm.is_open())
        throw std::invalid_argument("failed to open inputfile");
    
    std::string line;
    std::regex pattern(R"((\w+)\s*:\s*((\S+)))"); // var : data # ignore
    std::smatch matches;                          //            # comments
    while(std::getline(istrm, line)) // extract input file data
    {
        if(std::regex_search(line, matches, pattern))
        {
            std::string key = matches[1];
            std::string value = 
                matches[2];
            if(key == "rho_re") rho_re_file = value;
            else if(key == "rho_im") rho_im_file = value;
            else if(key == "hbar") hbar = std::stod(value);
            else if(key == "cap_gamma") cap_gamma = std::stod(value);
            else if(key == "cap_omega_plus") cap_omega_plus = std::stod(value);
            else if(key == "cap_omega_pi") cap_omega_pi = std::stod(value);
            else if(key == "cap_omega_minus") 
                cap_omega_minus = std::stod(value);
            else if(key == "cap_delta_B") cap_delta_B = std::stod(value);
            else if(key == "cap_delta_pi") cap_delta_pi = std::stod(value);
            else if(key == "cap_delta_plus") cap_delta_plus = std::stod(value);
            else if(key == "cap_delta_upper") 
                cap_delta_upper = std::stod(value);
            else if(key == "ti") ti = std::stod(value);
            else if(key == "tf") tf = std::stod(value);
        }
    }
    istrm.close();

    if(!rho_re_file.empty() && !rho_im_file.empty()) this->rho 
        = matrix(read(rho_re_file, ' ')) 
        + 1i*matrix(read(rho_im_file, ' '));

    // diagonal matrix elements
    H(0,0) = cap_delta_B - cap_delta_pi + cap_delta_plus;
    H(2,2) = -cap_delta_B - cap_delta_pi + cap_delta_plus;
    H(3,3) = cap_delta_upper - cap_delta_pi + cap_delta_plus; 

    // off-diagonal couplings
    H(0,3) = -cap_omega_plus/2; H(3,0) = -cap_omega_plus/2;
    H(1,3) = -cap_omega_pi/2; H(3,1) = -cap_omega_pi/2;
    H(2,3) = -cap_omega_minus/2; H(3,2) = -cap_omega_minus/2;
} // parameterized constructor from input file

// trying time dependent hamiltonian idea below
matrix<std::complex<double>> fourlevel_state::H2(double t)
{
double t_on = 3.5, t_off = 2.3, tau = 0.1;
double cap_omega_plus = 34.5575191895, cap_omega_pi = 15.7079632679,
       cap_omega_minus = 34.5575191895, cap_delta_B = 6.28318530718,
       cap_delta_pi = 12.5663706144, cap_delta_plus = 12.5663706144,
       cap_delta_upper = 0.;
    matrix<std::complex<double>> HH(4,4);
    
    // diagonal matrix elements
    HH(0,0) = cap_delta_B - cap_delta_pi + cap_delta_plus;
    HH(1,1) = 0.;
    HH(2,2) = -cap_delta_B - cap_delta_pi + cap_delta_plus;
    HH(3,3) = cap_delta_upper - cap_delta_pi + cap_delta_plus;

    // off-diagonal couplings. use if statements for now...
    if(t < t_off)
    {
        HH(0,3) = -cap_omega_plus/2; HH(3,0) = -cap_omega_plus/2;
        HH(2,3) = -cap_omega_minus/2; HH(3,2) = -cap_omega_minus/2;
    }
    else if(t >= t_off && t < (t_on + t_off)/2.)
    {
        HH(0,3) = -cap_omega_plus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
        HH(3,0) = -cap_omega_plus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
        HH(2,3) = -cap_omega_minus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
        HH(3,2) = -cap_omega_minus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
    }
    else if(t >= (t_on + t_off)/2. && t < t_on)
    {
        HH(0,3) = -cap_omega_plus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
        HH(3,0) = -cap_omega_plus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
        HH(2,3) = -cap_omega_minus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
        HH(3,2) = -cap_omega_minus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
    }
    else if(t >= t_on)
    {
        HH(0,3) = -cap_omega_plus/2; HH(3,0) = -cap_omega_plus/2;
        HH(2,3) = -cap_omega_minus/2; HH(3,2) = -cap_omega_minus/2;
    }

    HH(1,3) = -cap_omega_pi/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
//        - cap_omega_pi/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
    HH(3,1) = -cap_omega_pi/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
 //       - cap_omega_pi/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);

    return HH;
}

std::pair<std::vector<double>, 
std::vector<std::vector<std::complex<double>>>> fourlevel_state::solve()
{
    std::cerr << cap_gamma << '\n';
    auto master_eq = [&](double t, std::vector<std::complex<double>> rho_vector)
    {
        matrix<std::complex<double>> rho_matrix = stack(rho_vector, 4);
        matrix<std::complex<double>> sol
            = -1i*(H2(t)*rho_matrix - rho_matrix*H2(t));
/*        for(int i = 0; i < sol.size1; i++) sol(i,3) -= cap_gamma/2*rho(i,3);
        for(int i = 0; i < sol.size2; i++) sol(3,i) -= cap_gamma/2*rho(3,i);
        sol(0,0) += cap_gamma/3*rho(3,3);
        sol(1,1) += cap_gamma/3*rho(3,3);
        sol(2,2) += cap_gamma/3*rho(3,3);*/ // is this stupid?
        for(int i = 0; i < sol.size1; i++) 
            sol(i,3) -= cap_gamma/2*rho_matrix(i,3);
        for(int i = 0; i < sol.size2; i++) 
            sol(3,i) -= cap_gamma/2*rho_matrix(3,i);
        sol(0,0) += cap_gamma/3*rho_matrix(3,3);
        sol(1,1) += cap_gamma/3*rho_matrix(3,3);
        sol(2,2) += cap_gamma/3*rho_matrix(3,3);
        return sol.unravel();
    };
    return driver<std::complex<double>>(master_eq, {ti, tf}, rho.unravel());
}
/*
std::vector<cubic_spline<std::complex<double>>> fourlevel_state::solve()
{
    auto master_eq = [&](double t, std::vector<std::complex<double>> rho_vector)
    {
        matrix<std::complex<double>> rho_matrix = stack(rho_vector, 4);
        matrix<std::complex<double>> sol
            = -1i*(H*rho_matrix - rho_matrix*H);
        return sol.unravel();
    };
//    return driver<std::complex<double>>(master_eq, {ti, tf}, rho.unravel());
    return interpolant<std::complex<double>>(master_eq, {ti, tf}, rho.unravel());
}
*/
/*
std::vector<std::complex<double>> master_eq(double t, fourlevel_state state)
{
    matrix<std::complex<double>> sol 
        = -1i*(state.H*state.rho - state.rho*state.H);
    return sol.unravel();
} // master_eq*/
