#include <iostream>
#include "fourlevel.h"
#include "/home/asudler/git/capstone/misctools/misctools.h"

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

//    from_ifstream(inputfile);

    matrix<std::complex<double>> H(4,4), rho(4,4);
    double hbar = 3., cap_gamma = 4.;
    fourlevel_state myState(H, rho, hbar, cap_gamma);
    std::cout << myState.hbar << '\n';
    myState.rho.print();
//    myState.hbar = 1.61;
    std::cout << myState.hbar << '\n';

    fourlevel_state myState2(inputfile);
    std::cout << myState2.hbar << '\n';
    myState2.rho.print();

    return 0;
} // main

