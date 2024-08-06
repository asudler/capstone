#include <cmath>
#include <iostream>
#include "/home/ajs/repos/capstone/src/misctools/misctools.h"
#include "/home/ajs/repos/capstone/src/spline/spline.h"

void makedata(int argc, char *argv[])
{
    int ni = 0; // number of intervals in the data
    double xmax = 0.0; // max value for data

    for(int i = 0; i < argc; ++i)
    {
        if(std::string(argv[i]) == "-ni") ni = std::stoi(argv[i+1]);
        if(std::string(argv[i]) == "-xmax") xmax = std::stod(argv[i+1]);    
    }
    if(ni == 0 || xmax == 0.0)
        throw std::invalid_argument("makedata: invalid or missing "
            "cmd line inputs");

    for(int i = 0; i <= ni; i++)
    {
        double x = ((double)i/ni)*xmax;
        double y = std::exp(-0.1*x)*std::sin(x);
        std::cout << x << ',' << y << '\n';
    }
} // makedata

void test_cubic_spline(int argc, char *argv[])
{
    std::string fname = ""; // file name containing data to be splined
    int ni = 0; // number of intervals in the spline

    for(int i = 0; i < argc; ++i)
    {
        if(std::string(argv[i]) == "-fname") fname = argv[i+1];
        if(std::string(argv[i]) == "-ni") ni = std::stoi(argv[i+1]);
    }
    if(fname == "")
        throw std::invalid_argument("test_cubic_spline: "
            "input filename not properly provided");
    if(ni == 0)
        throw std::invalid_argument("test_cubic_spline: "
            "invalid or missing cmd line inputs");

    std::vector<std::vector<double>> v = csv_read::read(fname);
    std::vector<double> xs(v.size());
    std::vector<double> ys(v.size());
    for(int i = 0; i < v.size(); i++)
    {
        xs[i] = v[i][0]; ys[i] = v[i][1];
    } // extract x and y data from fileread
    
    cubic_spline cspline(xs, ys);

    std::cout << "x,f(x),f'(x),F(x)\n";
    for(int i = 0; i <= ni; i++)
    {
        double z = (xs[xs.size()-1]/ni)*i;
        std::cout << z << ',' << cspline.evaluate(z) << ',' << 
            cspline.derivative(z) << ',' << cspline.integral(z) << '\n';
    }
} // test_cubic_spline

int main(int argc, char *argv[]) 
{   
    for(int i = 0; i < argc; ++i)
    {
        if(std::string(argv[i]) == "--makedata") makedata(argc, argv);
        if(std::string(argv[i]) == "--test")
            test_cubic_spline(argc, argv);
    }
    return 0;
} // main

