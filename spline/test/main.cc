#include <cmath>
#include <complex>
#include <iostream>
#include "/home/asudler/git/capstone/misctools/misctools.h"
#include "/home/asudler/git/capstone/spline/spline.h"

using namespace std::complex_literals;

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

	std::cout << "x,Re(y),Im(y)\n";
    for(int i = 0; i <= ni; i++)
    {
        double x = ((double)i/ni)*xmax;
        std::complex<double> y = std::exp(-0.1*x)*std::sin(x)
			+ 1i*std::exp(-0.1*x)*std::cos(x);
        std::cout << x << ',' << y.real() << ',' << y.imag() << '\n';
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

    std::vector<std::vector<double>> v = read(fname);
    std::vector<double> xs(v.size());
    std::vector<std::complex<double>> ys(v.size());
    for(int i = 0; i < v.size(); i++)
    {
        xs[i] = v[i][0]; ys[i] = v[i][1] + 1i*v[i][2];
    } // extract x and y data from fileread
    
	cubic_spline<std::complex<double>> cspline(xs, ys);
    std::cout << "x,Re[f(x)],Im[f(x)],Re[f'(x)],Im[f'(x)],Re[F(x)],Im[F(x)]\n";
    for(int i = 0; i <= ni; i++)
    {
        double z = (xs[xs.size()-1]/ni)*i;
        std::cout << z << ',' << cspline.evaluate(z).real() << ',' <<
			cspline.evaluate(z).imag() << ',' << 
			cspline.derivative(z).real() << ',' <<
            cspline.derivative(z).imag() << ',' <<
			cspline.integral(z).real() << ',' << 
            cspline.integral(z).imag() << ',' << '\n';
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

