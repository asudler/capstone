#include <cmath>
#include <iostream>
#include <stdexcept>
#include "spline.h"

/* binsearch:
 * binary search algorithm; i.e.,
 * finding target interval of a sorted array
 * by bisection */
int spline::binsearch(const std::vector<double> x, double z)
{
    if(z < x[0] || z > x[x.size() - 1]) // reject searches outside array
        throw std::invalid_argument("binsearch: bad z");
    int i = 0; int j = x.size() - 1;
    while(j - i > 1)
    {
        int mid = (i+j)/2;
        if(z > x[mid]) i = mid; else j = mid;
    }
    return i; // lower or "left" index (beginning at 0)
} // binsearch

/* cubic_spline:
 * cubic spline interpolation class constructor */
cubic_spline::cubic_spline(const std::vector<double> xs,
    const std::vector<double> ys)
{
    if(xs.size() != ys.size())
        throw new std::invalid_argument("cubic_spline: x and y "
            "must be same dimension");
    x = xs; y = ys;
    
    int n = xs.size();
    b.resize(n); c.resize(n-1); d.resize(n-1);
    std::vector<double> h(n-1), p(n-1), D(n), Q(n-1), B(n);
    for(int i = 0; i < n-1; i++)
    {
        h[i] = x[i+1] - x[i];
        if(h[i] < 0) throw std::invalid_argument("cubic_spline: "
            "x must be sorted"); // x's must strictly increase
        p[i] = (y[i+1] - y[i])/h[i];
    }

    // building tridiagonal system
    D[0] = 2; Q[0] = 1; B[0] = 3*p[0];
    for(int i = 0; i < n - 2; i++) 
    {
        D[i+1] = 2*h[i]/h[i+1] + 2;
        Q[i+1] = h[i]/h[i+1];
        B[i+1] = 3*(p[i] + p[i+1]*h[i]/h[i+1]);
    }
    D[n-1] = 2; B[n-1] = 3*p[n-2];
    
    // gaussian elimination
    for(int i = 1; i < n; i++) 
    {
        D[i] -= Q[i-1]/D[i-1];
        B[i] -= B[i-1]/D[i-1];
    }
    
    // back substitution
    b[n-1] = B[n-1]/D[n-1];
    for(int i = n - 2; i >= 0; i--) 
        b[i] = (B[i] - Q[i]*b[i+1])/D[i];
    for(int i = 0; i < n - 1; i++) 
    {
        c[i] = (-2*b[i] - b[i+1] + 3*p[i])/h[i];
        d[i] = (b[i] + b[i+1] - 2*p[i])/(h[i]*h[i]);
    }
} // constructor

/* cubic_spline evaluate:
 * evaluate the cubic spline at a given coordinate input */
double cubic_spline::evaluate(double z) const
{
    int i = binsearch(x, z);
    return y[i] + b[i]*(z - x[i])
        + c[i]*std::pow(z - x[i], 2)
        + d[i]*std::pow(z - x[i], 3);
} // evaluate

/* cubic_spline derivative:
 * calculate the derivative of the cubic spline
 * at a given coordinate input */
double cubic_spline::derivative(double z) const
{
    int i = binsearch(x, z);
    return b[i] + 2*c[i]*(z - x[i])
        + 3*d[i]*std::pow(z - x[i], 2);
} // derivative

/* cubic_spline integral:
 * calculate the antiderivative of the cubic spline
 * from the beginning of the dataset up to
 * a given coordinate input */
double cubic_spline::integral(double z) const
{
    int i = binsearch(x, z);
    double integral = 0; double dx;
    for(int j = 0; j < i; j++) 
    {
        dx = x[j+1] - x[j];
        if(dx < 0) 
            throw std::invalid_argument("integral: x not sorted");
        integral += y[j]*dx + b[j]*std::pow(dx, 2)/2
            + c[j]*std::pow(dx, 3)/3 + d[j]*std::pow(dx, 4)/4;
    }
    integral += y[i]*(z - x[i]) + 1.0*b[i]*std::pow(z - x[i], 2)/2
        + c[i]*std::pow(z - x[i], 3)/3 + d[i]*std::pow(z - x[i], 4)/4;
    return integral;
} // integral

