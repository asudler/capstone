#ifndef SPLINE_H
#define SPLINE_H

#include <vector>

/* abstract base spline class 
 * could be extended to linear and quadratic splines in the future */
class spline
{
    public:
        spline() {}
        virtual ~spline() {}
    
    protected:
        static int binsearch(const std::vector<double>& x, double z);
}; // spline

/* cubic spline interpolation derived class
 * works with real and complex y data */
template <typename T>
class cubic_spline : public spline
{
    public:
        cubic_spline(const std::vector<double>& xs, const std::vector<T>& ys);
        virtual ~cubic_spline() {}
        virtual T evaluate(double z) const;
        virtual T derivative(double z) const;
        virtual T integral(double z) const;

    private:
        std::vector<double> x_;
        std::vector<T> y_, b_, c_, d_;
}; // cubic_spline

#endif // SPLINE_H

