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

// (1) initialize a state from input file parameters
// (2) solve the dynamics of this state
// (3) profit

matrix<std::complex<double>> fourlevel_state::rho_default_()
{
    matrix<std::complex<double>> result(4,4); // zeros
    return result;
} // default density matrix (4x4 of zeros)

fourlevel_state::fourlevel_state() : // default constructor
    rho(rho_default_()), hbar(1.), cap_gamma(1.),
    cap_omega_plus(0.), cap_omega_pi(0.), cap_omega_minus(0.), 
    cap_delta_B(0.), cap_delta_pi(0.), cap_delta_plus(0.),
    cap_delta_upper(0.), ti(0.), tf(1.), t_on(0.), t_off(0.), tau(0.) {}

fourlevel_state::fourlevel_state
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
) : // parameterized constructor
    rho(rho_init), hbar(hbar_init), cap_gamma(cap_gamma_init),
    cap_omega_plus(cap_omega_pi_init), cap_omega_pi(cap_omega_pi_init),
    cap_omega_minus(cap_omega_minus_init), cap_delta_B(cap_delta_B_init),
    cap_delta_pi(cap_delta_pi_init), cap_delta_plus(cap_delta_plus_init),
    cap_delta_upper(cap_delta_upper_init), ti(ti_init), tf(tf_init),
    t_on(t_on_init), t_off(t_off_init), tau(tau_init) {}

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
            std::string value = matches[2];
            if(key == "rho_re") 
                rho_re_file = value;
            else if(key == "rho_im") 
                rho_im_file = value;
            else if(key == "hbar") 
                hbar = std::stod(value);
            else if(key == "cap_gamma") 
                cap_gamma = std::stod(value);
            else if(key == "cap_omega_plus") 
                cap_omega_plus = std::stod(value);
            else if(key == "cap_omega_pi") 
                cap_omega_pi = std::stod(value);
            else if(key == "cap_omega_minus") 
                cap_omega_minus = std::stod(value);
            else if(key == "cap_delta_B") 
                cap_delta_B = std::stod(value);
            else if(key == "cap_delta_pi") 
                cap_delta_pi = std::stod(value);
            else if(key == "cap_delta_plus") 
                cap_delta_plus = std::stod(value);
            else if(key == "cap_delta_upper") 
                cap_delta_upper = std::stod(value);
            else if(key == "ti") 
                ti = std::stod(value);
            else if(key == "tf") 
                tf = std::stod(value);
            else if(key == "t_on")
                t_on = std::stod(value);
            else if(key == "t_off")
                t_off = std::stod(value);
            else if(key == "tau")
                tau = std::stod(value);
        }
    }
    istrm.close(); // perhaps there is a better way to do this

    if(!rho_re_file.empty() && !rho_im_file.empty())
        rho = matrix(read(rho_re_file, ' ')) 
            + 1i*matrix(read(rho_im_file, ' '));
} // parameterized constructor from input file

/* this gives the time-dependent Hamiltonian for the system 
 * note: the current process could be somewhat inefficient 
 * since it needlessly recaluates time-independent values each time.
 * REFACTOR?  */
matrix<std::complex<double>> fourlevel_state::H(double t)
{
    matrix<std::complex<double>> M(4,4);
    
    // diagonal matrix elements
    // not time dependent
    M(0,0) = cap_delta_B - cap_delta_pi + cap_delta_plus;
    M(1,1) = 0.;
    M(2,2) = -cap_delta_B - cap_delta_pi + cap_delta_plus;
    M(3,3) = cap_delta_upper - cap_delta_pi + cap_delta_plus;

    // off-diagonal couplings. use if statements for now...
    if(t < t_on)
    {
        M(0,3) = -cap_omega_plus/2; M(3,0) = -cap_omega_plus/2;
        M(2,3) = -cap_omega_minus/2; M(3,2) = -cap_omega_minus/2;
    }
    else if(t >= t_on && t < (t_off + t_on)/2.)
    {
        M(0,3) = -cap_omega_plus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
        M(3,0) = -cap_omega_plus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
        M(2,3) = -cap_omega_minus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
        M(3,2) = -cap_omega_minus/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
    }
    else if(t >= (t_off + t_on)/2. && t < t_off)
    {
        M(0,3) = -cap_omega_plus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
        M(3,0) = -cap_omega_plus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
        M(2,3) = -cap_omega_minus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
        M(3,2) = -cap_omega_minus/2*exp(-0.5*(t - t_off)*(t - t_off)/tau/tau);
    }
    else if(t >= t_off)
    {
        M(0,3) = -cap_omega_plus/2; M(3,0) = -cap_omega_plus/2;
        M(2,3) = -cap_omega_minus/2; M(3,2) = -cap_omega_minus/2;
    }

    M(1,3) = -cap_omega_pi/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);
    M(3,1) = -cap_omega_pi/2*exp(-0.5*(t - t_on)*(t - t_on)/tau/tau);

    return M;
} // H(t) (time-dependent Hamiltonian)

std::pair<std::vector<double>, 
std::vector<std::vector<std::complex<double>>>> fourlevel_state::solve()
{
    auto master_eq = [&](double t, std::vector<std::complex<double>> rho_vector)
    {
        matrix<std::complex<double>> rho_matrix = stack(rho_vector, 4);
        matrix<std::complex<double>> sol
            = -1i*(H(t)*rho_matrix - rho_matrix*H(t));
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
} // equation of motion for state
