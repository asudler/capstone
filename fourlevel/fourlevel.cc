#include <cmath>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include "fourlevel.h"
#include "/home/asudler/git/capstone/linalg/matrix/matrix.h"
#include "/home/asudler/git/capstone/misctools/misctools.h"
#include "/home/asudler/git/capstone/rkf45/rkf45.h"
#include "/home/asudler/git/capstone/spline/spline.h"

using namespace std::complex_literals;

// (1) initialize a state from input file parameters
// (2) solve the dynamics of this state
// (3) profit

filenames::filenames() : beams(""), rho(""), rho_log(""), spatial_log(""), 
    spatial_rho_xii_log(""), spatial_rho_xif_log(""), spatial_omega_base(""),
    physical_grid(""), polariton_grid(""), physical_grid_check(""), 
    physical_grid_lite(""), omega_print(0) {}


filenames::filenames(std::string inputfile) :
filenames()
{
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
            if(key == "beams") 
                beams = value;
            else if(key == "rho") 
                rho = value;
            else if(key == "rho_log") 
                rho_log = value;
            else if(key == "spatial_log") 
                spatial_log = value;
            else if(key == "spatial_rho_xii_log") 
                spatial_rho_xii_log = value;
            else if(key == "spatial_rho_xif_log")
                spatial_rho_xif_log = value;
            else if(key == "spatial_omega_base") 
                spatial_omega_base = value;
            else if(key == "omega_print") 
                omega_print = std::stoi(value);
            else if(key == "polaritons_base")
                polaritons_base = value;
            else if(key == "physical_grid")
                physical_grid = value;
            else if(key == "physical_grid_lite")
                physical_grid_lite = value;
            else if(key == "polariton_grid")
                polariton_grid = value;
            else if(key == "physical_grid_check")
                physical_grid_check = value;
        }
    }
    istrm.close();
}


matrix<std::complex<double>> fourlevel_state::rho0_default_()
{
    return matrix<std::complex<double>>(4,4); // zeros
} // default density matrix (4x4 of zeros)


cubic_spline<std::complex<double>> fourlevel_state::cap_omega_spline_default_()
{
    std::vector<double> ts = {0., 1.}; // no default ctr for cubic_spline :(
    std::vector<std::complex<double>> ys = {0., 0.};
    return cubic_spline<std::complex<double>>(ts, ys);
} // default spline for cap_omega_pi_t


fourlevel_state::fourlevel_state() : // default constructor
    rho0(rho0_default_()), hbar(1.), cap_gamma(1.),
    cap_omega_plus(0.), cap_omega_pi(0.), cap_omega_minus(0.), 
    cap_delta_B(0.), cap_delta_pi(0.), cap_delta_plus(0.),
    cap_delta_upper(0.), cap_omega_plus_t(cap_omega_spline_default_()), 
    cap_omega_pi_t(cap_omega_spline_default_()),
    cap_omega_minus_t(cap_omega_spline_default_()) {}


fourlevel_state::fourlevel_state
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
) : // parameterized constructor
    rho0(rho0_init), hbar(hbar_init), cap_gamma(cap_gamma_init),
    cap_omega_plus(cap_omega_pi_init), cap_omega_pi(cap_omega_pi_init),
    cap_omega_minus(cap_omega_minus_init), cap_delta_B(cap_delta_B_init),
    cap_delta_pi(cap_delta_pi_init), cap_delta_plus(cap_delta_plus_init),
    cap_delta_upper(cap_delta_upper_init), 
    cap_omega_plus_t(cap_omega_minus_t_init),
    cap_omega_pi_t(cap_omega_pi_t_init),
    cap_omega_minus_t(cap_omega_plus_t_init) {}


fourlevel_state::fourlevel_state(std::string inputfile) : fourlevel_state()
{
    std::string rho0_re_file, rho0_im_file;
    
    std::ifstream istrm(inputfile.c_str());
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
            if(key == "rho0_re") 
                rho0_re_file = value;
            else if(key == "rho0_im") 
                rho0_im_file = value;
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
            else if(key == "decoherence") 
                decoherence = std::stod(value);
            else if(key == "ti") 
                ti = std::stod(value);
            else if(key == "tf") 
                tf = std::stod(value);
            else if(key == "dt") 
                dt = std::stod(value);
            else if(key == "nt") 
                nt = std::stoi(value);
            else if(key == "const_dt")
                std::istringstream(value) >> std::boolalpha >> const_dt;
            else if(key == "t_on_pi") 
                t_on_pi = std::stod(value);
            else if(key == "t_off_pi") 
                t_off_pi = std::stod(value);
            else if(key == "tau_pi") 
                tau_pi = std::stod(value);
            else if(key == "t_on1_pm") 
                t_on1_pm = std::stod(value);
            else if(key == "t_off1_pm") 
                t_off1_pm = std::stod(value);  
            else if(key == "t_on2_pm") 
                t_on2_pm = std::stod(value);
            else if(key == "t_off2_pm") 
                t_off2_pm = std::stod(value);
            else if(key == "tau_pm") 
                tau_pm = std::stod(value);
            else if(key == "g")
                g = std::stod(value);
            else if(key == "nn")
                nn = std::stod(value);
            else if(key == "chi_m")
                chi_m = std::stod(value);
            else if(key == "chi_p")
                chi_p = std::stod(value);
            else if(key == "use_tanh")
                use_tanh = std::stoi(value);
            else if(key == "rise1")
                rise1 = std::stod(value);
            else if(key == "fall1")
                fall1 = std::stod(value);
            else if(key == "rise2")
                rise2 = std::stod(value);
            else if(key == "fall2")
                fall2 = std::stod(value);
            else if(key == "t_phase")
                t_phase = std::stod(value);
        }
    }
    istrm.close(); // perhaps there is a better way to do this

    if(!rho0_re_file.empty() && !rho0_im_file.empty())
        rho0 = matrix(read(rho0_re_file, ' ')) 
            + 1i*matrix(read(rho0_im_file, ' '));

    matrix<std::complex<double>> decoherence_help(4,4);
    decoherence_help(0,1) = decoherence;
    decoherence_help(0,2) = decoherence;
    decoherence_help(0,3) = decoherence;
    decoherence_help(1,0) = decoherence;
    decoherence_help(1,2) = decoherence;
    decoherence_help(1,3) = decoherence;
    decoherence_help(2,0) = decoherence;
    decoherence_help(2,1) = decoherence;
    decoherence_help(2,3) = decoherence;
    decoherence_help(3,0) = decoherence;
    decoherence_help(3,1) = decoherence;
    decoherence_help(3,2) = decoherence;
    decoherence_matrix = decoherence_help;

    // create control beam phase functions
//    f_chi_p = [&](double t)
//    {
//        if(t < t_phase) return 0.;
//        else if(t >= t_phase) return chi_p;
//        else return 0.;
//    };
//    f_chi_m = [&](double t)
//    {
//        if(t < t_phase) return 0.;
//        else if(t >= t_phase) return chi_m;
//        else return 0.;
//    };

    // initialize beam profiles
    std::vector<double> ts; double t = ti;
    while(t < tf) { ts.push_back(t); t += 0.01; }
    ts.push_back(tf);

    std::vector<std::complex<double>> spline_help_plus, spline_help_pi,
        spline_help_minus;
    for(double t : ts) 
    {
        std::complex<double> help = cap_omega_pi
            *exp(-0.5*(t - t_on_pi)*(t - t_on_pi)/tau_pi/tau_pi);
        spline_help_pi.push_back(help);

        if(use_tanh == 1)
        {
            help = 0.5*cap_omega_plus*(
                    std::tanh(rise1*(t - t_on1_pm))
                    - std::tanh(fall1*(t - t_off1_pm))
                    + (std::tanh(rise2*(t - t_on2_pm))
                    - std::tanh(fall2*(t - t_off2_pm)))*std::exp(1i*chi_p)
                   );
//            help *= std::exp(1i*f_chi_p(t));
            spline_help_plus.push_back(help);
            help = 0.5*cap_omega_minus*(
                    std::tanh(rise1*(t - t_on1_pm))
                    - std::tanh(fall1*(t - t_off1_pm))
                    + (std::tanh(rise2*(t - t_on2_pm))
                    - std::tanh(fall2*(t - t_off2_pm)))*std::exp(1i*chi_m)
                   );
//            help *= std::exp(1i*f_chi_m(t));
            spline_help_minus.push_back(help);
        }
        else 
        {
            if(t < t_on1_pm)
            {
                help = cap_omega_plus
                    *exp(-0.5*(t - t_on1_pm)*(t - t_on1_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_p(t));
                spline_help_plus.push_back(help);
                help = cap_omega_minus
                    *exp(-0.5*(t - t_on1_pm)*(t - t_on1_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_m(t));
                spline_help_minus.push_back(help);
            }
            else if(t >= t_on1_pm && t < t_off1_pm)
            {
                spline_help_plus.push_back(cap_omega_plus
                    *std::exp(1i*f_chi_p(t)));
                spline_help_minus.push_back(cap_omega_minus
                    *std::exp(1i*f_chi_m(t)));
            }
            else if(t >= t_off1_pm && t < (t_off1_pm + t_on2_pm)/2.)
            {
                help = cap_omega_plus
                    *exp(-0.5*(t - t_off1_pm)*(t - t_off1_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_p(t));
                spline_help_plus.push_back(help);
                help = cap_omega_minus
                    *exp(-0.5*(t - t_off1_pm)*(t - t_off1_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_m(t));
                spline_help_minus.push_back(help);
            }
            else if(t >= (t_off1_pm + t_on2_pm)/2. && t < t_on2_pm)
            {
                help = cap_omega_plus
                    *exp(-0.5*(t - t_on2_pm)*(t - t_on2_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_p(t));
                spline_help_plus.push_back(help);
                help = cap_omega_minus
                    *exp(-0.5*(t - t_on2_pm)*(t - t_on2_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_m(t));
                spline_help_minus.push_back(help);
            }
            else if(t >= t_on2_pm && t < t_off2_pm)
            {
                spline_help_plus.push_back(cap_omega_plus
                    *std::exp(1i*f_chi_p(t)));
                spline_help_minus.push_back(cap_omega_minus
                    *std::exp(1i*f_chi_m(t)));
            }
            else if(t >= t_off2_pm)
            {
                help = cap_omega_plus
                    *exp(-0.5*(t - t_off2_pm)*(t - t_off2_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_p(t));
                spline_help_plus.push_back(help);
                help = cap_omega_minus
                    *exp(-0.5*(t - t_off2_pm)*(t - t_off2_pm)
                    /tau_pm/tau_pm);
                help *= std::exp(1i*f_chi_m(t));
                spline_help_minus.push_back(help);
            }
        }
    }

    cap_omega_plus_t 
        = cubic_spline<std::complex<double>>(ts, spline_help_plus);
    cap_omega_pi_t 
        = cubic_spline<std::complex<double>>(ts, spline_help_pi);
    cap_omega_minus_t 
        = cubic_spline<std::complex<double>>(ts, spline_help_minus);

    // create mixing angle functions from beam splines
    /*
    theta = [&](double t) 
    { 
        return std::atan(std::sqrt(g*g*nn/
               (std::abs(cap_omega_plus_t.evaluate(t))
                *std::abs(cap_omega_plus_t.evaluate(t))
               + std::abs(cap_omega_minus_t.evaluate(t))
                *std::abs(cap_omega_minus_t.evaluate(t)))));
    };*/
    theta = [&](double t) 
    { 
        return std::atan2(
            std::sqrt(g*g*nn),
            std::sqrt(std::abs(cap_omega_plus_t.evaluate(t))
            *std::abs(cap_omega_plus_t.evaluate(t))
            + std::abs(cap_omega_minus_t.evaluate(t))
            *std::abs(cap_omega_minus_t.evaluate(t)))
        );
    };
    phi = [&](double t)
    {
        return std::atan2(std::abs(cap_omega_minus_t.evaluate(t)),
               std::abs(cap_omega_plus_t.evaluate(t)));
    };
} // parameterized constructor from input file


/* boundary conditions struct initialization below 
 * it might be better to put this in separate .cc and .h files */
boundary_conditions::boundary_conditions() : // default ctr
    mu_alpha(0.), epsilon(0.7), tolerance(1.) {}

boundary_conditions::boundary_conditions(std::string inputfile) :
boundary_conditions()
{
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
            if(key == "xi_min") 
                xi_min = std::stod(value);
            else if(key == "xi_max") 
                xi_max = std::stod(value);
            else if(key == "nxi") 
                nxi = std::stoi(value);
            else if(key == "mu_alpha") 
                mu_alpha = std::stod(value);
            else if(key == "epsilon") 
                epsilon = std::stod(value);
            else if(key == "tolerance") 
                tolerance = std::stod(value);
            else if(key == "nn")
                nn = std::stoi(value);
        }
    }   
    istrm.close(); // once again maybe there's a more efficient way here...
}


/* this gives the time-dependent Hamiltonian for the system 
 * note: the current process could be somewhat inefficient 
 * since it needlessly recaluates time-independent values each time.
 * REFACTOR? */
matrix<std::complex<double>> fourlevel_state::H(double t)
{
    matrix<std::complex<double>> M(4,4);
    
    // diagonal matrix elements
    // not time dependent
    M(0,0) = cap_delta_B - cap_delta_pi + cap_delta_plus;
    M(1,1) = 0.;
    M(2,2) = -cap_delta_B - cap_delta_pi + cap_delta_plus;
    M(3,3) = cap_delta_upper - cap_delta_pi + cap_delta_plus;

    // off-diagonal couplings
    M(0,3) = -cap_omega_plus_t.evaluate(t)/2.; M(3,0) = conj(M(0,3));
    M(2,3) = -cap_omega_minus_t.evaluate(t)/2.; M(3,2) = conj(M(2,3));
        
    // coupling between |2> and auxillary state is z-dependent
    M(1,3) = -cap_omega_pi_t.evaluate(t)/2.; M(3,1) = conj(M(1,3));

    return M;
} // H(t) (time-dependent Hamiltonian)


/* rotate:
 * matrix for rotating from physical to polariton picture */
matrix<std::complex<double>> fourlevel_state::rotate(double t)
{
    matrix<std::complex<double>> M(3,3);

    M(0,0) = std::cos(theta(t))/g/std::sqrt(nn);
    M(0,1) = -std::exp(1i*f_chi_m(t))*std::sin(theta(t))*std::sin(phi(t));
    M(0,2) = -std::exp(1i*f_chi_p(t))*std::sin(theta(t))*std::cos(phi(t));
    M(1,0) = std::sin(theta(t))/g/std::sqrt(nn);
    M(1,1) = std::exp(1i*f_chi_m(t))*std::cos(theta(t))*std::sin(phi(t));
    M(1,2) = std::exp(1i*f_chi_p(t))*std::cos(theta(t))*std::cos(phi(t));
    M(2,0) = 0;
    M(2,1) = std::exp(-1i*f_chi_p(t))*std::cos(phi(t));
    M(2,2) = -std::exp(-1i*f_chi_m(t))*std::sin(phi(t));

    return M;
}


/* unrotate:
 * matrix for rotating from polariton to physical picture */
matrix<std::complex<double>> fourlevel_state::unrotate(double t)
{
    matrix<std::complex<double>> M(3,3);

    M(0,0) = g*std::sqrt(nn)*std::cos(theta(t));
    M(0,1) = g*std::sqrt(nn)*std::sin(theta(t));
    M(0,2) = 0;
    M(1,0) = -std::exp(-1i*f_chi_m(t))*std::sin(theta(t))*std::sin(phi(t));
    M(1,1) = std::exp(-1i*f_chi_m(t))*std::cos(theta(t))*std::sin(phi(t));
    M(1,2) = std::exp(1i*f_chi_p(t))*std::cos(phi(t));
    M(2,0) = -std::exp(-1i*f_chi_p(t))*std::sin(theta(t))*std::cos(phi(t));
    M(2,1) = std::exp(-1i*f_chi_p(t))*std::cos(theta(t))*std::cos(phi(t));
    M(2,2) = -std::exp(1i*f_chi_m(t))*std::cos(phi(t));

    return M;
}


// ARCHAIC AND WRONG! DO NOT USE
void fourlevel_state::print_beams(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "BEAM PROFILES\n"
      << "t, plus beam (re), plus beam (im), "
      << "pi beam (re), pi beam (im), "
      << "minus beam (re), minus beam (im)\n";

    double t = ti;
    while(t < tf)
    {
        f << t << ' '
          << cap_omega_plus_t.evaluate(t).real() << ' '
          << cap_omega_plus_t.evaluate(t).imag() << ' '
          << cap_omega_pi_t.evaluate(t).real() << ' '
          << cap_omega_pi_t.evaluate(t).imag() << ' '
          << cap_omega_minus_t.evaluate(t).real() << ' '
          << cap_omega_minus_t.evaluate(t).imag() << '\n';
        t += 0.01;
    } // WHAT ARE WE PRINTING? INTENSITY? --> should be squared
      // RABI COUPLINGS? --> should not be called "plus beam (re)" eg
      // BE EXACT!!!

    f.close();
}

void fourlevel_state::print_rabi_couplings(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "Rabi couplings.\n"
      << "t, Re(\\Omega_+), Im(\\Omega_+), "
      << "Re(\\Omega_{\\pi}), Im(\\Omega_{\\pi}), "
      << "Re(\\Omega_-), Im(\\Omega_-), "
      << "\\theta, \\phi\n";

    double t = ti;
    while(t < tf)
    {
        f << t << ' '
          << cap_omega_plus_t.evaluate(t).real() << ' '
          << cap_omega_plus_t.evaluate(t).imag() << ' '
          << cap_omega_pi_t.evaluate(t).real() << ' '
          << cap_omega_pi_t.evaluate(t).imag() << ' '
          << cap_omega_minus_t.evaluate(t).real() << ' '
          << cap_omega_minus_t.evaluate(t).imag() << ' '
          << theta(t) << ' ' << phi(t) << '\n';
        t += 0.002;
    }

    f.close();
}


std::vector<std::complex<double>> fourlevel_state::master_eq_
(
    double t,
    std::vector<std::complex<double>> rho_vector
)
{
    matrix<std::complex<double>> rho_matrix = stack(rho_vector, 4);
    matrix<std::complex<double>> sol
        = -1i*(H(t)*rho_matrix - rho_matrix*H(t));
    for(int i = 0; i < sol.size1; i++)
    {
        for(int j = 0; j < sol.size2; j++)
        {
            sol(i,j) -= decoherence_matrix(i,j)*rho_matrix(i,j);
        }
    }
    for(int i = 0; i < sol.size1; i++)
        sol(i,3) -= cap_gamma/2.*rho_matrix(i,3);
    for(int i = 0; i < sol.size2; i++)
        sol(3,i) -= cap_gamma/2.*rho_matrix(3,i);
    if(cap_omega_minus != 0. && cap_omega_plus != 0.)
    {
        sol(0,0) += cap_gamma/3.*rho_matrix(3,3);
        sol(1,1) += cap_gamma/3.*rho_matrix(3,3);
        sol(2,2) += cap_gamma/3.*rho_matrix(3,3);
    }
    else if(cap_omega_minus == 0. && cap_omega_plus != 0.) 
    {
        sol(0,0) += cap_gamma/2.*rho_matrix(3,3);
        sol(1,1) += cap_gamma/2.*rho_matrix(3,3);
    }
    else if(cap_omega_minus != 0. && cap_omega_plus == 0.)
    {
        sol(1,1) += cap_gamma/2.*rho_matrix(3,3);
        sol(2,2) += cap_gamma/2.*rho_matrix(3,3);
    }
    else
    {
        sol(1,1) += cap_gamma*rho_matrix(3,3);
    }
    return sol.unravel();
} // master_eq, i.e. the equation of motion for the density matrix rho


void fourlevel_state::solve()
{
    auto rhs = [&](double t, std::vector<std::complex<double>> rho_vector)
    { return master_eq_(t, rho_vector); };

    if(dt == 0 && nt == 0)
    {
        if(const_dt) throw std::invalid_argument("const_dt cannot be set to "
            "true if both dt and nt are set to 0");
        auto [t, rho]
           = driver<std::complex<double>>(rhs, {ti, tf}, rho0.unravel());
        times = t;
        solutions = rho;
    }
    else 
    {   
        if(nt != 0) dt = tf/(nt - 1.); // reassign dt if nonzero nt
        if(const_dt) // use simple RK4 routine
        {
            auto [t, rho] = driver<std::complex<double>>(rhs, {ti, tf}, 
                rho0.unravel(), dt, true);
            times = t;
            solutions = rho;
        }
        else // use RKF45 and return the splined solution eval'ed at times
        {
            auto sol = interpolant<std::complex<double>>(rhs, {ti, tf}, 
                rho0.unravel(), 0.125, false, 
                std::numeric_limits<double>::infinity(), 1e-9, 1e-9);

            // get output
            std::vector<std::vector<std::complex<double>>> ys;
            std::vector<double> ts; 
            double t = ti;
            while(t < tf)
            {
                if(t + dt > tf) t = tf;
                ts.push_back(t);

                std::vector<std::complex<double>> help;
                for(int i = 0; i < sol.size(); i++)
                {
                    help.push_back(sol[i].evaluate(t));
                }
                ys.push_back(help);

                t += dt;
            }

            times = ts;
            solutions = ys;
        }
    }

//    // also spline the rho solution
//    for(std::vector<std::complex<double>> sol : solutions)
//    {
//        auto help_spline = cubic_spline<std::complex<double>>(times, sol);
//        solutions_spline.push_back(help
//    }
} // solve, i.e. solution via numerical ode methods


void fourlevel_state::print_rho(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "DENSITY MATRIX VS TIME\n"
      << "t, rho_{ij}.real(), rho_{ij}.imag()\n";

    for(int i = 0; i < times.size(); i++)
    {
        f << times[i] << ' ';
        for(int j = 0; j < solutions[i].size(); j++)
        {
            f << solutions[i][j].real() << ' ' 
              << solutions[i][j].imag() << ' ';
        }
        f << '\n';
    }

    f.close();
}


void fourlevel_state::print_rho_log(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "DENSITY MATRIX VS TIME\n"
      << "time, [matrix elem.], rho (re), rho (im)\n\n";

    for(int i = 0; i < times.size(); i++)
    {
        int row = 1, column = 1;
        for(int j = 0; j < solutions[i].size(); j++)
        {
            f << times[i] << '\t'
              << "[" << row << column << "]" << ' '
              << solutions[i][j].real() << ' ' 
              << solutions[i][j].imag() << '\n';
            column += 1;
            if((j+1) % rho0.size1 == 0) { row += 1; column = 1; }
        }
        f << '\n';
    }

    f.close();
}


void fourlevel_state::print_polaritons(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "PSI, PHI, and Z polaritons as functions of TIME\n"
      << "t, Re(Psi), Im(Psi), Re(Phi), Im(Phi), Re(Z), Im(Z)\n";
    
    for(int i = 0; i < times.size(); i++)
    {
        f << times[i] << ' ';

        // make a vector containing E (Omega_{pi}/g), rho23, and rho21
        std::vector<std::complex<double>> physical_vector
            = {cap_omega_pi_t.evaluate(times[i])/g, solutions[i][5],
               solutions[i][7]};

        // make a new vector containing Psi, Phi, and Z
        auto polariton_vector = rotate(times[i])*physical_vector;

        // print to file
        for(int j = 0; j < polariton_vector.size(); j++)
        {
            f << polariton_vector[j].real() << ' ' 
              << polariton_vector[j].imag() << ' ';
        }
        f << '\n';
    }

    f.close();
}


void fourlevel_state::print_polaritons_log(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "POLARITONS vs TIME\n"
      << "time, [polariton], Re(polariton), Im(polariton)\n";
    
    for(int i = 0; i < times.size(); i++)
    {
        f << times[i] << '\t';

        // make a vector containing E (Omega_{pi}/g), rho23, and rho21
        std::vector<std::complex<double>> physical_vector
            = {cap_omega_pi_t.evaluate(times[i])/g, solutions[i][5],
               solutions[i][7]};

        // make a new vector containing Psi, Phi, and Z
        auto polariton_vector = rotate(times[i])*physical_vector;

        // print to file
        f << "[Psi] " 
          << polariton_vector[0].real() << ' ' 
          << polariton_vector[0].imag() << '\n'
          << "[Phi] " 
          << polariton_vector[1].real() << ' ' 
          << polariton_vector[1].imag() << '\n'
          << "[Z] " 
          << polariton_vector[2].real() << ' ' 
          << polariton_vector[2].imag() << '\n';
    }

    f.close();
}
