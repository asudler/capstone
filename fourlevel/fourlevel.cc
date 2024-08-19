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

matrix<std::complex<double>> fourlevel_state::H_default_()
{
    matrix<std::complex<double>> result(4,4);
    result(3,3) = 1.; // excited state energy = 1
    return result;
} // default Hamiltonian

matrix<std::complex<double>> fourlevel_state::rho_default_()
{
    matrix<std::complex<double>> result(4,4);
    result(0,0) = 1.;
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
            else if(key == "ti") ti = std::stod(value);
            else if(key == "tf") tf = std::stod(value);
        }
    }
    istrm.close();

    if(!rho_re_file.empty() && !rho_im_file.empty()) this->rho 
        = matrix(read(rho_re_file, ' ')) 
        + 1i*matrix(read(rho_im_file, ' '));
} // parameterized constructor from input file

std::pair<std::vector<double>, 
std::vector<std::vector<std::complex<double>>>> fourlevel_state::solve()
{
    auto master_eq = [&](double t, std::vector<std::complex<double>> rho_vector)
    {
        matrix<std::complex<double>> rho_matrix = stack(rho_vector, 4);
        matrix<std::complex<double>> sol
            = -1i*(H*rho_matrix - rho_matrix*H);
        return sol.unravel();
    };
    return driver<std::complex<double>>(master_eq, {ti, tf}, rho.unravel());
}
/*
std::vector<std::complex<double>> master_eq(double t, fourlevel_state state)
{
    matrix<std::complex<double>> sol 
        = -1i*(state.H*state.rho - state.rho*state.H);
    return sol.unravel();
} // master_eq*/
