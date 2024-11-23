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
    omega_print(0) {}


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
        throw std::invalid_argument("FFfailed to open inputfile");
    
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
        }
    }
    istrm.close(); // perhaps there is a better way to do this

    if(!rho0_re_file.empty() && !rho0_im_file.empty())
        rho0 = matrix(read(rho0_re_file, ' ')) 
            + 1i*matrix(read(rho0_im_file, ' '));

    // initialize beam profiles
    std::vector<double> ts; double t = 0;
    while(t < tf) { ts.push_back(t); t += 0.01; }
    ts.push_back(tf);

    std::vector<std::complex<double>> spline_help_plus, spline_help_pi,
        spline_help_minus;
    for(double t : ts) 
    {
        std::complex<double> help = cap_omega_pi
            *exp(-0.5*(t - t_on_pi)*(t - t_on_pi)/tau_pi/tau_pi);
        spline_help_pi.push_back(help);

        if(t < t_on1_pm)
        {
            help = cap_omega_plus
                *exp(-0.5*(t - t_on1_pm)*(t - t_on1_pm)
                /tau_pm/tau_pm);
            spline_help_plus.push_back(help);
            help = cap_omega_minus
                *exp(-0.5*(t - t_on1_pm)*(t - t_on1_pm)
                /tau_pm/tau_pm);
            spline_help_minus.push_back(help);
        }
        else if(t >= t_on1_pm && t < t_off1_pm)
        {
            spline_help_plus.push_back(cap_omega_plus);
            spline_help_minus.push_back(cap_omega_minus);
        }
        else if(t >= t_off1_pm && t < (t_off1_pm + t_on2_pm)/2.)
        {
            help = cap_omega_plus
                *exp(-0.5*(t - t_off1_pm)*(t - t_off1_pm)
                /tau_pm/tau_pm);
            spline_help_plus.push_back(help);
            help = cap_omega_minus
                *exp(-0.5*(t - t_off1_pm)*(t - t_off1_pm)
                /tau_pm/tau_pm);
            spline_help_minus.push_back(help);
        }
        else if(t >= (t_off1_pm + t_on2_pm)/2. && t < t_on2_pm)
        {
            help = cap_omega_plus
                *exp(-0.5*(t - t_on2_pm)*(t - t_on2_pm)
                /tau_pm/tau_pm);
            spline_help_plus.push_back(help);
            help = cap_omega_minus
                *exp(-0.5*(t - t_on2_pm)*(t - t_on2_pm)
                /tau_pm/tau_pm);
            spline_help_minus.push_back(help);
        }
        else if(t >= t_on2_pm && t < t_off2_pm)
        {
            spline_help_plus.push_back(cap_omega_plus);
            spline_help_minus.push_back(cap_omega_minus);
        }
        else if(t >= t_off2_pm)
        {
            help = cap_omega_plus
                *exp(-0.5*(t - t_off2_pm)*(t - t_off2_pm)
                /tau_pm/tau_pm);
            spline_help_plus.push_back(help);
            help = cap_omega_minus
                *exp(-0.5*(t - t_off2_pm)*(t - t_off2_pm)
                /tau_pm/tau_pm);
            spline_help_minus.push_back(help);
        }
    }

    cap_omega_plus_t 
        = cubic_spline<std::complex<double>>(ts, spline_help_plus);
    cap_omega_pi_t 
        = cubic_spline<std::complex<double>>(ts, spline_help_pi);
    cap_omega_minus_t 
        = cubic_spline<std::complex<double>>(ts, spline_help_minus);
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


void fourlevel_state::print_beams(std::string file)
{
    std::ofstream f;
    f.open(file);
    f << "BEAM PROFILES\n"
      << "t, plus beam (re), plus beam (im), "
      << "pi beam (re), pi beam (im), "
      << "minus beam (re), minus beam (im)\n";

    double t = 0;
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
        sol(i,3) -= cap_gamma/2.*rho_matrix(i,3);
    for(int i = 0; i < sol.size2; i++)
        sol(3,i) -= cap_gamma/2.*rho_matrix(3,i);
    sol(0,0) += cap_gamma/3.*rho_matrix(3,3);
    sol(1,1) += cap_gamma/3.*rho_matrix(3,3);
    sol(2,2) += cap_gamma/3.*rho_matrix(3,3);
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

