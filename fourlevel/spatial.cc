#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "fourlevel.h"
#include "spatial.h"
#include "/home/asudler/git/capstone/linalg/matrix/matrix.h"
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
        // read \Omega_{\pi} from file
        // needs to be tested
        // if(true /*change later*/)
        // {
        //     auto 
        // } 
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
        */
        if(files.spatial_rho_xif_log != "")
            states0[bc.nxi - 1].print_rho_log(files.spatial_rho_xif_log);
        if(norm < bc.tolerance) // consider "normalizing" lol
        { 
            return states0; 
        }

        omega0 = omega1.copy();
        states0 = states1;
    }
} // spatial_solve
