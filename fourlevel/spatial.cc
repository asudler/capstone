#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "/home/asudler/git/capstone/fourlevel/fourlevel.h"
#include "spatial.h"
#include "/home/asudler/git/capstone/linalg/matrix/matrix.h"
#include "/home/asudler/git/capstone/misctools/misctools.h"
#include "/home/asudler/git/capstone/spline/spline.h"

using namespace std::complex_literals;

std::vector<fourlevel_state> spatial_solve
(
    fourlevel_state state,
    boundary_conditions bc, // for the z-grid
    filenames files // for monitoring
)
{
    std::ofstream f;
    if(files.spatial_log != "")
    {
        f.open(files.spatial_log);
        f.close();
    }

    if(files.spatial_omega_base != "") // if writing omega to file,
    {                                  // make sure omega_print != 0
        if(files.omega_print == 0)     // otherwise divide by 0 error
            throw std::invalid_argument("spatial_omega_base exists, so "
                "ensure omega_print is not set to 0");
    }
    
    std::vector<fourlevel_state> states0(bc.nxi);
    for(int i = 0; i < bc.nxi; i++)
    {
        states0[i] = state;        
        states0[i].solve();
    }

    if(files.spatial_rho_xii_log != "") 
        states0[0].print_rho_log(files.spatial_rho_xii_log);

    matrix<std::complex<double>> omega0(bc.nxi, states0[0].times.size());
    std::vector<double> ts = states0[0].times;
    for(int i = 0; i < bc.nxi; i++)
    {
        for(int j = 0; j < ts.size(); j++)
        {
            omega0(i,j) = states0[i].cap_omega_pi_t.evaluate(ts[j]);
        }
    }

    // sum the magnitude of rho at all grid points for the 0th iteration
    double metric_re = 0, metric_im = 0;
    for(int i = 0; i < bc.nxi; i++)
    {
        for(int j = 0; j < ts.size(); j++)
        {
            for(int k = 0; k < states0[i].rho0.size1*states0[i].rho0.size2; k++)
            {
                metric_re += std::abs(states0[i].solutions[j][k].real());
                metric_im += std::abs(states0[i].solutions[j][k].imag());
            }
        }
    }

    if(files.spatial_log != "") // print magnitude of rho across the grid
    {
        f.open(files.spatial_log, std::ios_base::app); // append
        f << std::setprecision(16)
          << "0" << '\t' << "METRIC_RE: " << metric_re << '\n'
          << "0" << '\t' << "METRIC_IM: " << metric_im << '\n'
          << '\n';
        f.close();
    }
    /* actually, we don't want to do this beyond testing purposes (probably)
    if(files.spatial_omega_base != "") // write omega0 to file
    {
        std::ofstream omega_file_re;
        omega_file_re.open(files.spatial_omega_base + "_0000_re.dat");
        std::ofstream omega_file_im;
        omega_file_im.open(files.spatial_omega_base + "_0000_im.dat");
        real(omega0).print("pi-beam spatial matrix (re)", ' ', 
            omega_file_re);
        imag(omega0).print("pi-beam spatial matrix (im)", ' ', 
            omega_file_im);
        omega_file_re.close();
        omega_file_im.close();
    }
    */

    double dxi = (bc.xi_max - bc.xi_min)/(bc.nxi - 1.);
    int counter = 0;
    while(true)
    {
        counter++;
        
        auto omega1 = omega0.copy();
        for(int i = 0; i < bc.nxi - 1; i++)
        {
            for(int j = 0; j < ts.size(); j++)
            {
                omega1(i+1,j) = omega0(i,j)
                    - 1i*dxi*bc.mu_alpha*states0[i].solutions[j][7]; // rho24
                //  ^ "+" sign used in Shan calculations 2025-02-04.
            }
        }

        omega1 = bc.epsilon*omega0 + (1. - bc.epsilon)*omega1;
    
        auto states1 = states0;
        double norm = 0, norm_re = 0, norm_im = 0;
        metric_re = 0, metric_im = 0;
        for(int i = 0; i < bc.nxi; i++)
        {
            states1[i].cap_omega_pi_t 
                = cubic_spline<std::complex<double>>(ts, omega1.row(i));
            states1[i].solve();

            // measure convergence
            for(int j = 0; j < ts.size(); j++)
            {
                for(int k = 0; 
                    k < states0[i].rho0.size1*states0[i].rho0.size2; 
                    k++)
                {
                    metric_re += std::abs(states1[i].solutions[j][k].real());
                    metric_im += std::abs(states1[i].solutions[j][k].imag());
                    norm += std::abs(states1[i].solutions[j][k]
                         - states0[i].solutions[j][k]);
                    norm_re += std::abs((states1[i].solutions[j][k]
                         - states0[i].solutions[j][k]).real());
                    norm_im += std::abs((states1[i].solutions[j][k]
                         - states0[i].solutions[j][k]).imag());
                }
            }
        }

        if(files.spatial_log != "")
        {
            f.open(files.spatial_log, std::ios_base::app); // append
            f << std::setprecision(16)
              << counter << '\t' << "METRIC_RE: " << metric_re << "\n"
              << counter << '\t' << "METRIC_IM: " << metric_im << "\n"
              << counter << '\t' << "DIFF_RE: " << norm_re << "\n"
              << counter << '\t' << "DIFF_IM: " << norm_im << "\n"
              << counter << '\t' << "DIFF: " << norm << "\n"
              << "\n";    
            f.close();
        }
        /* again, probably don't want this...
        if(counter % files.omega_print == 0) // print omega1 to files
        {                                    // both real and imag parts
            std::string index = std::to_string(counter);
            while(index.length() < 4) index.insert(0, "0");

            std::ofstream omega_file_re; 
            omega_file_re.open(files.spatial_omega_base + "_" + 
                  index + "_re.dat");
            std::ofstream omega_file_im; 
            omega_file_im.open(files.spatial_omega_base + "_" + 
                  index + "_im.dat");
            real(omega1).print("pi-beam spatial matrix (re)", ' ', 
                omega_file_re);
            imag(omega1).print("pi-beam spatial matrix (im)", ' ', 
                omega_file_im);
            omega_file_re.close();
            omega_file_im.close();
        }
        // end long comment here*/
        if(files.spatial_rho_xif_log != "")
            states0[bc.nxi - 1].print_rho_log(files.spatial_rho_xif_log);
        if(norm < bc.tolerance) // consider "normalizing" lol
        { 
            return states0; 
        }

        omega0 = omega1.copy();
        states0 = states1;
    }
    // for testing purposes
    //return states0;
} // spatial_solve

std::vector<fourlevel_state> spatial_solve
(
    std::vector<fourlevel_state> state,
    boundary_conditions bc, // for the z-grid
    filenames files // for monitoring
)
{
    std::ofstream f;
    if(files.spatial_log != "")
    {
        f.open(files.spatial_log);
        f.close();
    }
    
    if(files.spatial_omega_base != "") // if writing omega to file,
    {                                  // make sure omega_print != 0
        if(files.omega_print == 0)     // otherwise divide by 0 error
            throw std::invalid_argument("spatial_omega_base exists, so "
                "ensure omega_print is not set to 0");
    }
    
    std::vector<fourlevel_state> states0(bc.nxi);
    std::cerr << "state size: " << state.size() << '\n';
    std::cerr << "states0 size: " << states0.size() << '\n';
    for(int i = 0; i < bc.nxi; i++)
    {
        std::cerr << i << ' ';
        states0[i] = state[i]; // use converged state from func args
        states0[i].solve();
    }

    if(files.spatial_rho_xii_log != "") 
        states0[0].print_rho_log(files.spatial_rho_xii_log);

    matrix<std::complex<double>> omega0(bc.nxi, states0[0].times.size());
    std::vector<double> ts = states0[0].times;
    for(int i = 0; i < bc.nxi; i++)
    {
        for(int j = 0; j < ts.size(); j++)
        {
            omega0(i,j) = states0[i].cap_omega_pi_t.evaluate(ts[j]);
        }
    }

    // sum the magnitude of rho at all grid points for the 0th iteration
    double metric_re = 0, metric_im = 0;
    for(int i = 0; i < bc.nxi; i++)
    {
        for(int j = 0; j < ts.size(); j++)
        {
            for(int k = 0; k < states0[i].rho0.size1*states0[i].rho0.size2; k++)
            {
                metric_re += std::abs(states0[i].solutions[j][k].real());
                metric_im += std::abs(states0[i].solutions[j][k].imag());
            }
        }
    }

    if(files.spatial_log != "") // print magnitude of rho across the grid
    {
        f.open(files.spatial_log, std::ios_base::app); // append
        f << std::setprecision(16)
          << "0" << '\t' << "METRIC_RE: " << metric_re << '\n'
          << "0" << '\t' << "METRIC_IM: " << metric_im << '\n'
          << '\n';
        f.close();
    }

    double dxi = (bc.xi_max - bc.xi_min)/(bc.nxi - 1.);
    int counter = 0;
    while(true)
    {
        counter++;
        
        auto omega1 = omega0.copy();
        for(int i = 0; i < bc.nxi - 1; i++)
        {
            for(int j = 0; j < ts.size(); j++)
            {
                omega1(i+1,j) = omega0(i,j)
                    - 1i*dxi*bc.mu_alpha*states0[i].solutions[j][7]; // rho24
                //  ^ "+" sign used in Shan calculations 2025-02-04.
            }
        }

        omega1 = bc.epsilon*omega0 + (1. - bc.epsilon)*omega1;
    
        auto states1 = states0;
        double norm = 0, norm_re = 0, norm_im = 0;
        metric_re = 0, metric_im = 0;
        for(int i = 0; i < bc.nxi; i++)
        {
            states1[i].cap_omega_pi_t 
                = cubic_spline<std::complex<double>>(ts, omega1.row(i));
            states1[i].solve();

            // measure convergence
            for(int j = 0; j < ts.size(); j++)
            {
                for(int k = 0; 
                    k < states0[i].rho0.size1*states0[i].rho0.size2; 
                    k++)
                {
                    metric_re += std::abs(states1[i].solutions[j][k].real());
                    metric_im += std::abs(states1[i].solutions[j][k].imag());
                    norm += std::abs(states1[i].solutions[j][k]
                         - states0[i].solutions[j][k]);
                    norm_re += std::abs((states1[i].solutions[j][k]
                         - states0[i].solutions[j][k]).real());
                    norm_im += std::abs((states1[i].solutions[j][k]
                         - states0[i].solutions[j][k]).imag());
                }
            }
        }

        if(files.spatial_log != "")
        {
            f.open(files.spatial_log, std::ios_base::app); // append
            f << std::setprecision(16)
              << counter << '\t' << "METRIC_RE: " << metric_re << "\n"
              << counter << '\t' << "METRIC_IM: " << metric_im << "\n"
              << counter << '\t' << "DIFF_RE: " << norm_re << "\n"
              << counter << '\t' << "DIFF_IM: " << norm_im << "\n"
              << counter << '\t' << "DIFF: " << norm << "\n"
              << "\n";    
            f.close();
        }

        if(files.spatial_rho_xif_log != "")
            states0[bc.nxi - 1].print_rho_log(files.spatial_rho_xif_log);
        if(norm < bc.tolerance)
        { 
            return states0; 
        }

        omega0 = omega1.copy();
        states0 = states1;
    }
    // for testing purposes
    //return states0;
} // spatial_solve


/* spatial_print:
 * print a nice gnuplot friendly grid
 * of the physical picture parameters */
void spatial_print
(
    std::vector<fourlevel_state> states,
    boundary_conditions bc,
    filenames files
)
{
    // Format:
    // \xi', \tau, Re(\Omega_{\pi}), Im(\Omega_{\pi}),
    // Re(\rho_{ij}), Im(rho_{ij}) (i,j=1..4)
    std::ofstream f;
    f.precision(15);
    f.open(files.physical_grid);
    f << "\\xi', \\tau, Re(\\Omega_{\\pi}), Im(\\Omega_{\\pi}), "
      << "Re(\\rho_{ij}), Im(rho_{ij}) (i,j=1..4)\n\n";
    
    double dxi = (bc.xi_max - bc.xi_min)/(bc.nxi - 1.);
    double xi = 0.;
    
    for(auto state : states)
    {
        for(int ti = 0; ti < state.times.size(); ti++)
        {
            f << xi << '\t' << state.times[ti] << '\t'
              << state.cap_omega_pi_t.evaluate(state.times[ti]).real() << '\t'
              << state.cap_omega_pi_t.evaluate(state.times[ti]).imag() << '\t';
            for(int i = 0; i < state.rho0.size1*state.rho0.size2; i++)
            {
                f << state.solutions[ti][i].real() << '\t'
                  << state.solutions[ti][i].imag() << '\t';
            }
            f << '\n';
        }
        f << '\n';
        xi += dxi;
    }

    // job done! Hopefully! :)
    f.close();
}


/* spatial_print_lite:
 * print a nice gnuplot friendly grid
 * of the physical picture parameters 
 * but also just some coherences of rho */
void spatial_print_lite
(
    std::vector<fourlevel_state> states,
    boundary_conditions bc,
    filenames files
)
{
    // Format:
    // \xi', \tau, Re(\Omega_{\pi}), Im(\Omega_{\pi}),
    // Re(\rho_{ij}), Im(rho_{ij}) (just 21,23,13)
    std::ofstream f;
    f.precision(15);
    f.open(files.physical_grid_lite);
    f << "\\xi', \\tau, Re(\\Omega_{\\pi}), Im(\\Omega_{\\pi}), "
      << "Re(\\rho_{}), Im(rho_{ij}) (ij=21,23,13)\n\n";
    
    double dxi = (bc.xi_max - bc.xi_min)/(bc.nxi - 1.);
    double xi = 0.;
    
    for(auto state : states)
    {
        for(int ti = 0; ti < state.times.size(); ti++)
        {
            f << xi << '\t' << state.times[ti] << '\t'
              << state.cap_omega_pi_t.evaluate(state.times[ti]).real() << '\t'
              << state.cap_omega_pi_t.evaluate(state.times[ti]).imag() << '\t'
              << state.solutions[ti][5].real() << '\t'
              << state.solutions[ti][5].imag() << '\t'
              << state.solutions[ti][7].real() << '\t'
              << state.solutions[ti][7].imag() << '\t'
              << state.solutions[ti][3].real() << '\t'
              << state.solutions[ti][3].imag() << '\t'
              << '\n';
        }
        f << '\n';
        xi += dxi;
    }

    // job done! Hopefully! :)
    f.close();
}


/* spatial_read:
 * generate a grid of spatial states
 * from a file formatted the way that the spatial_print function
 * outputs its data */
std::vector<fourlevel_state> spatial_read
(
    std::string inputfile,
    boundary_conditions bc,
    filenames files
)
{
    auto data = read(files.physical_grid, '\t', 2);
    double tf = data[data.size() - 1][1];
    std::vector<double> times_help;
    std::vector<std::complex<double>> cap_omega_pi_help;
    std::vector<std::vector<std::complex<double>>> solutions_help;
    std::vector<fourlevel_state> states;

    for(auto line : data)
    {
        times_help.push_back(line[1]);
        cap_omega_pi_help.push_back(line[2] + 1i*line[3]);
        std::vector<std::complex<double>> solutions_row_help;
        for(int i = 0; i < 16; i++)
        {
            solutions_row_help.push_back(line[4 + 2*i] + 1i*line[5 + 2*i]);
        }
        solutions_help.push_back(solutions_row_help);

        if(line[1] == tf)
        {
            fourlevel_state state(inputfile);
            state.times = times_help;
            state.solutions = solutions_help;
            state.cap_omega_pi_t = cubic_spline<std::complex<double>>
                (state.times, cap_omega_pi_help);
            states.push_back(state);
            
            // Remove items from vectors
            times_help.clear();
            solutions_help.clear();
            cap_omega_pi_help.clear();
        }
    }
    
    return states;
}


/* rotate:
 * read spatially-dependent physical data
 * rotate to polariton picture
 * print result */
void rotate(fourlevel_state state, filenames files)
{
    auto physical_zt = read(files.physical_grid, '\t', 1);
    std::ofstream f;
    f.open(files.polariton_grid);
    f << "t,z,Re[\\Psi/sqrt(nn)],Im[\\Psi/sqrt(nn)],"
      << "Re[\\Phi/sqrt(nn)],Im[\\Psi/sqrt(nn)],"
      << "Re[Z/sqrt(nn)],Im[Z/sqrt(nn)]\n";

    double z = 0;
    for(auto line : physical_zt)
    {
        if(line[0] != z) f << '\n';
        std::vector<std::complex<double>> rhs = {
            line[2] + 1i*line[3], line[4] + 1i*line[5], line[6] + 1i*line[7]
        };

        z = line[0];
        auto lhs = state.rotate(line[1])*rhs;
        f << z << '\t' << line[1] << '\t'
          << lhs[0].real() << '\t' << lhs[0].imag() << '\t'
          << lhs[1].real() << '\t' << lhs[1].imag() << '\t'
          << lhs[2].real() << '\t' << lhs[2].imag() << '\n';
    }
    f.close();
}


/* unrotate:
 * read spatially-dependent polariton data
 * rotate to physical picture
 * print result */
void unrotate(fourlevel_state state, filenames files)
{
    auto polariton_zt = read(files.polariton_grid, '\t', 1);
    std::ofstream f;
    f.open(files.physical_grid_check);
    f << "t,z,Re[\\Omega_{\\pi}],Im[\\Omega_{\\pi}],"
      << "Re[\\rho_{23}],Im[\\rho_{23}],"
      << "Re[\\rho_{21}],Im[\\rho_{21}]\n";

    double z = 0;
    for(auto line : polariton_zt)
    {
        if(line[0] != z) f << '\n';
        std::vector<std::complex<double>> rhs = {
            line[2] + 1i*line[3], line[4] + 1i*line[5], line[6] + 1i*line[7]
        };
        
        z = line[0];
        auto lhs = state.unrotate(line[1])*rhs;
        f << z << '\t' << line[1] << '\t'
          << lhs[0].real() << '\t' << lhs[0].imag() << '\t'
          << lhs[1].real() << '\t' << lhs[1].imag() << '\t'
          << lhs[2].real() << '\t' << lhs[2].imag() << '\n';
    }
    f.close();
}
