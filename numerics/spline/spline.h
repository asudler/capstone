#ifndef SPLINE_H
#define SPLINE_H

#include <vector>

// abstract base spline class
class spline
{
    public:
        spline() {}
        virtual ~spline() {}

        // virtual methods for derived classes
        virtual double evaluate(double z) const = 0;
        virtual double derivative(double z) const = 0;
        virtual double integral(double z) const = 0;
    
    protected:
        static int binsearch(const std::vector<double> x, double z);
}; // spline

// cubic spline interpolation derived class
class cubic_spline : public spline
{
    public:
        cubic_spline(const std::vector<double> xs,
            const std::vector<double> ys);
        virtual ~cubic_spline() {}
        virtual double evaluate(double z) const override;
        virtual double derivative(double z) const override;
        virtual double integral(double z) const override;

    private:
        std::vector<double> x, y, b, c, d;
};

#endif // SPLINE_H

