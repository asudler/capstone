#include <algorithm>
#include <cmath>
#include <complex>

// perhaps the norm functions should be in different files
// since they are not specific to runge-kutta routines

/* norm:
 * calculate the vector norm of a vector of doubles */
double norm(std::vector<double> vec)
{
    double meanabs = 0;
    for(int i = 0; i < vec.size(); i++) meanabs += std::abs(vec[i]);
    if(meanabs == 0) meanabs = 1;
    meanabs /= vec.size();
    double sum = 0;
    for(int i = 0; i < vec.size(); i++)
        sum += (vec[i]/meanabs)*(vec[i]/meanabs);
    return meanabs*std::sqrt(sum);
} // norm

/* norm:
 * calculate the vector norm of a vector of complex doubles */
double norm(std::vector<std::complex<double>> vec)
{
    double meanabs = 0;
    for(int i = 0; i < vec.size(); i++) meanabs += std::abs(vec[i]);
    if(meanabs == 0) meanabs = 1;
    meanabs /= vec.size();
    double sum = 0;
    for(int i = 0; i < vec.size(); i++)
        sum += std::norm(vec[i])/meanabs/meanabs;
    return meanabs*std::sqrt(sum);
} // norm


/* rkf45:
 * adaptive runge-kutta 4th-order method with error estimate
 * https://en.wikipedia.org/wiki/Runge–Kutta–Fehlberg_method */
template <typename T>
std::pair<std::vector<T>, std::vector<T>> rkf45
(
    std::function<std::vector<T>(double, std::vector<T>)> f,
    double x,
    std::vector<T> y,
    double h
)
{
    int n = y.size();
    std::vector<double> a0 = {1./4.};
    std::vector<double> a1 = {3./32., 9./32.};
    std::vector<double> a2 = {1932./2197., -7200./2197., 7296./2197.};
    std::vector<double> a3 = {439./216., -8., 3680./513., -845./4104.};
    std::vector<double> a4
        = {-8./27., 2., -3544./2565., 1859./4104., -11./40.};
    std::vector<double> b0
        = {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.};
    std::vector<double> b1
        = {25./216., 0., 1408./2565., 2197./4104., -1./5., 0.};
    std::vector<double> c
        = {0., 1./4., 3./8., 12./13., 1., 1./2.};
    std::vector<T> yt(n), yh(n), yhs(n), yerr(n);

    /* question for the future:
     * would it be faster if all of these were in the same loop? */

    // k0 step
    std::vector<T> k0 = f(x,y);
    for(int i = 0; i < n; i++) yt[i] = y[i] + a0[0]*k0[i]*h;

    // k1 step
    std::vector<T> k1 = f(x + c[1]*h, yt);
    for(int i = 0; i < n; i++) yt[i] = y[i] + (a1[0]*k0[i] + a1[1]*k1[i])*h;

    // k2 step
    std::vector<T> k2 = f(x + c[2]*h, yt);
    for(int i = 0; i < n; i++) yt[i] = y[i] + (a2[0]*k0[i] + a2[1]*k1[i]
        + a2[2]*k2[i])*h;

    // k3 step
    std::vector<T> k3 = f(x + c[3]*h, yt);
    for(int i = 0; i < n; i++) yt[i] = y[i] + (a3[0]*k0[i] + a3[1]*k1[i]
        + a3[2]*k2[i] + a3[3]*k3[i])*h;

    // k4 step
    std::vector<T> k4 = f(x + c[4]*h, yt);
    for(int i = 0; i < n; i++) yt[i] = y[i] + (a4[0]*k0[i] + a4[1]*k1[i]
        + a4[2]*k2[i] + a4[3]*k3[i] + a4[4]*k4[i])*h;

    // k5 step and error
    std::vector<T> k5 = f(x + c[5]*h, yt);
    for(int i = 0; i < n; i++)
    {
        yh[i] = y[i] + (b0[0]*k0[i] + b0[2]*k2[i] + b0[3]*k3[i]
            + b0[4]*k4[i] + b0[5]*k5[i])*h;
        yhs[i] = y[i] + (b1[0]*k0[i] + b1[2]*k2[i] + b1[3]*k3[i]
            + b1[4]*k4[i] + b1[5]*k5[i])*h;
        yerr[i] = yh[i] - yhs[i];
    }

    return {yh, yerr};
} // rkf45


/* rk4:
 * fourth-order runga-kutta method
 * nothing fancy */
template <typename T>
std::vector<T> rk4
(
    std::function<std::vector<T>(double, std::vector<T>)> f,
    double x,
    std::vector<T> y,
    double h
)
{
    int n = y.size();
    std::vector<T> yt(n), yh(n);

    std::vector<T> k1 = f(x,y); // k1 step
    for(int i = 0; i < n; i++) { k1[i] *= h; yt[i] = y[i] + 0.5*k1[i]; }

    std::vector<T> k2 = f(x + 0.5*h, yt); // k2 step
    for(int i = 0; i < n; i++) { k2[i] *= h; yt[i] = y[i] + 0.5*k2[i]; }

    std::vector<T> k3 = f(x + 0.5*h, yt); // k3 step
    for(int i = 0; i < n; i++) { k3[i] *= h; yt[i] = y[i] + k3[i]; }

    std::vector<T> k4 = f(x + h, yt); // k4 step
    for(int i = 0; i < n; i++)
    {
        k4[i] *= h;
        yh[i] = y[i] + (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i])/6.;
    }

    return yh;
}


/* driver:
 * solve ODE with an adaptive step size 
 * uses rkf45 */
template <typename T>
std::pair<std::vector<double>, std::vector<std::vector<T>>> driver
(
    std::function<std::vector<T>(double, std::vector<T>)> f,
    std::pair<double, double> interval,
    std::vector<T> yi,
    double h,
    bool const_h,
    double max_h,
    double acc,
    double eps
)
{
    auto [a, b] = interval;
    double x = a; std::vector<T> y = yi;
    std::vector<double> xs; xs.push_back(x);
    std::vector<std::vector<T>> ys; ys.push_back(y);

    if(!const_h)
    {
        do // main driver loop
        {
            if(x >= b) return {xs, ys}; // job done, break loop
            if(x + h > b) h = b - x; // last step should end at b
            auto [yh, yerr] = rkf45(f, x, y, h);
            double tol = (acc + eps*norm(yh))*std::sqrt(h/(b-a));
            double err = norm(yerr);
            if(err <= tol) // accept step
            {
                x += h; y = yh;
                xs.push_back(x); ys.push_back(y);
            }
            h *= std::min(std::pow(tol/err, 0.25)*0.95, 2.); // adjust step
            if(h > max_h) h = max_h;
        }
        while(true);
    }
    else // constant step --> use simple rk4 with step = h
    {
        do // here max_h does not even matter, the only step size will be h
        {
            if(x + h > b) return {xs, ys}; // job done, break loop
            auto yh = rk4(f, x, y, h);
            x += h; y = yh;
            xs.push_back(x); ys.push_back(y);
        }
        while(true);
    }
} // driver


/* interpolant:
 * returns cubic spline interpolant of driver data */
template <typename T>
std::vector<cubic_spline<T>> interpolant
(
    const std::function<std::vector<T>(double, std::vector<T>)> f,
    const std::pair<double, double> interval,
    const std::vector<T> yi,
    double h,
    bool const_h,
    double max_h,
    double acc,
    double eps
)
{
    auto [x, ys] = driver(f, interval, yi, h, const_h, max_h, acc, eps);
    /* in the future,
     * consider allocating memory for the output
     * instead of pushing back (could be faster) */
    std::vector<cubic_spline<T>> output;
    for(int i = 0; i < ys[0].size(); i++) // interpolate all y solutions
    {
        std::vector<T> y(x.size());
        for(int j = 0; j < ys.size(); j++) y[j] = ys[j][i];
        cubic_spline<T> cspline(x, y);
        output.push_back(cspline);
    }
    return output;
} // interpolant
