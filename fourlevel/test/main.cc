#include <iostream>
#include "/home/asudler/git/capstone/fourlevel/fourlevel.h"
#include "/home/asudler/git/capstone/misctools/misctools.h"
#include "/home/asudler/git/capstone/rkf45/rkf45.h"

int main(int argc, char* argv[])
{
    std::string inputfile = ""; // input filename
    for(int i = 0; i < argc; ++i)
    {
        if(std::string(argv[i]) == "--inputfile" 
            || std::string(argv[i]) == "-i") inputfile = argv[i+1];
    }
    if(inputfile == "") throw std::invalid_argument("main: "
        "Input file not provided. Use --inputfile <file> or -i <file> "
        "when executing the program.");
/*
    matrix<std::complex<double>> H(4,4), rho(4,4);
    double hbar = 3., cap_gamma = 4.;
    fourlevel_state myState(H, rho, hbar, cap_gamma);
    std::cout << "hbar " << myState.hbar << '\n';
    myState.rho.print();

    fourlevel_state myState2(inputfile);
    std::cout << "hbar " << myState2.hbar << '\n';
    myState2.rho.print();

    fourlevel_state myState3(inputfile);
    std::cout << "hbar " << myState3.hbar << '\n';
    myState3.rho.print();

    fourlevel_state state;

    std::pair<std::vector<double>, 
        std::vector<std::vector<std::complex<double>>>> sol
        = driver<std::complex<double>>(master_eq, {state.ti, state.tf},
          state.rho); */
    
    fourlevel_state state;

    matrix<std::complex<double>> rho0(4,4);
    rho0(3,3) = 1;
    state.rho = rho0;

    matrix<std::complex<double>> H0(4,4);
    H0(0,0) = -0.1; H0(2,2) = 0.1; H0(3,3) = 1;
    H0(0,3) = -0.5; H0(3,0) = -0.5;
    H0(1,3) = -0.4; H0(3,1) = -0.5;
    H0(2,3) = -0.3; H0(3,2) = -0.3;
    state.H = H0;
    
    state.tf = 10.;

    auto [t, rho] = state.solve();
    for(int i = 0; i < rho.size(); i++)
    {
        std::cout << t[i] << ' ';
        for(int j = 0; j < rho[0].size(); j++)
        {
            std::cout << rho[i][j].real() << ' ' << rho[i][j].imag() << ' ';
        }
        std::cout << '\n';
    }
    return 0;
} // main

