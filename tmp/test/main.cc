#include <iostream>
#include "/home/asudler/git/capstone/tmp/fourlevel.h"
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
        "when executing the program."); // annoying: refactor into try/catch
   
    fourlevel_state state(inputfile);
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

