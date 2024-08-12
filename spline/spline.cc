#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include "spline.h"

/* explicit instantiation of template instances */
template class cubic_spline<double>;
template class cubic_spline<std::complex<double>>;

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
template <typename T>
cubic_spline<T>::cubic_spline(const std::vector<double> &xs, 
	const std::vector<T> &ys)
{
    if(xs.size() != ys.size())
        throw new std::invalid_argument("cubic_spline: x's and y's "
            "must have the same size");
    x_ = xs; y_ = ys;

    int n = x_.size();
    b_.resize(n); c_.resize(n-1); d_.resize(n-1);
    std::vector<double> h(n-1);
	std::vector<T> p(n-1), D(n), Q(n-1), B(n);
    for(int i = 0; i < n-1; i++)
    {
        h[i] = x_[i+1] - x_[i];
        if(h[i] < 0) throw std::invalid_argument("cubic_spline: "
            "x must be sorted"); // x's must strictly increase
        p[i] = (y_[i+1] - y_[i])/h[i];
    }

    // building tridiagonal system
    D[0] = 2.; Q[0] = 1.; B[0] = 3.*p[0];
    for(int i = 0; i < n - 2; i++) 
    {
        D[i+1] = 2.*h[i]/h[i+1] + 2.;
        Q[i+1] = h[i]/h[i+1];
        B[i+1] = 3.*(p[i] + p[i+1]*h[i]/h[i+1]);
    }
    D[n-1] = 2.; B[n-1] = 3.*p[n-2];

    // gaussian elimination
    for(int i = 1; i < n; i++) 
    {
        D[i] -= Q[i-1]/D[i-1];
        B[i] -= B[i-1]/D[i-1];
    }
    
    // back substitution
    b_[n-1] = B[n-1]/D[n-1];
    for(int i = n - 2; i >= 0; i--) 
        b_[i] = (B[i] - Q[i]*b_[i+1])/D[i];
    for(int i = 0; i < n - 1; i++) 
    {
        c_[i] = (-2.*b_[i] - b_[i+1] + 3.*p[i])/h[i];
        d_[i] = (b_[i] + b_[i+1] - 2.*p[i])/(h[i]*h[i]);
    }
} // constructor

/* cubic_spline evaluate:
 * evaluate the cubic spline at a given coordinate input */
template <typename T>
T cubic_spline<T>::evaluate(double z) const
{
    int i = binsearch(x_, z);
    return y_[i] + b_[i]*(z - x_[i])
        + c_[i]*std::pow(z - x_[i], 2)
        + d_[i]*std::pow(z - x_[i], 3);
} // evaluate

/* cubic_spline derivative:
 * calculate the derivative of the cubic spline
 * at a given coordinate input */
template <typename T>
T cubic_spline<T>::derivative(double z) const
{
    int i = binsearch(x_, z);
    return b_[i] + 2.*c_[i]*(z - x_[i]) + 3.*d_[i]*std::pow(z - x_[i], 2);
} // derivative

/* cubic_spline integral:
 * calculate the antiderivative of the cubic spline
 * from the beginning of the dataset up to
 * a given coordinate input */
template <typename T>
T cubic_spline<T>::integral(double z) const
{
    int i = binsearch(x_, z);
    T integral = 0; double dx;
    for(int j = 0; j < i; j++) 
    {
        dx = x_[j+1] - x_[j];
        integral += y_[j]*dx + b_[j]*std::pow(dx, 2)/2.
            + c_[j]*std::pow(dx, 3)/3. 
			+ d_[j]*std::pow(dx, 4)/4.;
    }
    integral += y_[i]*(z - x_[i]) 
		+ b_[i]*std::pow(z - x_[i], 2)/2. 
		+ c_[i]*std::pow(z - x_[i], 3)/3. 
		+ d_[i]*std::pow(z - x_[i], 4)/4.;
    return integral;
} // integral

