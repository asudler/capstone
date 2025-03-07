#include <iostream>
#include "/home/asudler/git/capstone/fourlevel/fourlevel.h"
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
    boundary_conditions bc(inputfile);
    filenames files(inputfile);

    // check to make sure the beam sequence is correct
    state.print_beams(files.beams);

    // get a solution
    state.solve();

    // print the solution
    state.print_rho(files.rho);

    // print an easier-to-debug log of the solution
    state.print_rho_log(files.rho_log);

    return 0;
} // main

