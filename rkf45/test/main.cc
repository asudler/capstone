#include <cmath>
#include <iostream>
#include "/home/asudler/git/capstone/src/rkf45/rkf45.h"

std::vector<double> f(double x, std::vector<double> y)
{
    return {y[1], -0.25*y[1] - 5.0*std::sin(y[0])};
} // f

void test_rkf45(int argc, char *argv[])
{   
    double b = 0, c = 0, yi = 0, dyi = 0, ti = 0, tf = 0,
        dt = 0; // cmd line inputs for oscillator system
    for(int i = 0; i < argc; ++i)
    {
        if(std::string(argv[i]) == "-b") b = std::stod(argv[i+1]);
        if(std::string(argv[i]) == "-c") c = std::stod(argv[i+1]);
        if(std::string(argv[i]) == "-yi") yi = std::stod(argv[i+1]);
        if(std::string(argv[i]) == "-dyi") dyi = std::stod(argv[i+1]);
        if(std::string(argv[i]) == "-ti") ti = std::stod(argv[i+1]);
        if(std::string(argv[i]) == "-tf") tf = std::stod(argv[i+1]);
        if(std::string(argv[i]) == "-dt") dt = std::stod(argv[i+1]);
    }

    std::vector<double> y0 = {yi, dyi};
    std::cout << "t,y,dy\n";

    if(dt == 0)
    {
        std::pair<std::vector<double>, std::vector<std::vector<double>>> sol
            = driver(f, {ti, tf}, y0);
        for(int i = 0; i < sol.first.size(); i++)
            std::cout << sol.first[i] << ',' << sol.second[i][0] << ',' 
                << sol.second[i][1] << '\n';
    }
    else
    {
        std::vector<cubic_spline> ys = interpolant(f, {ti, tf}, y0);
        while(ti < tf)
        {
            std::cout << ti << ',';
            for(int i = 0; i < ys.size(); i++)
                std::cout << ys[i].evaluate(ti) << ',';
            std::cout << '\n';
            ti += dt;
        }
        std::cout << tf << ',';
        for(int i = 0; i < ys.size(); i++) 
            std::cout << ys[i].evaluate(tf) << ',';
        std::cout << '\n';
    }
} // test_rkf45

int main(int argc, char *argv[])
{   
    for(int i = 0; i < argc; ++i)
    {
        if(std::string(argv[i]) == "--test")
            test_rkf45(argc, argv);
    }
    return 0;
} // main
