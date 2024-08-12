#ifndef RKF45_H
#define RKF45_H

#include <functional>
#include <utility>
#include <vector>
#include "/home/asudler/git/capstone/src/spline/spline.h"

std::pair<std::vector<double>, std::vector<double>> rkf45
(
    const std::function<std::vector<double>(double, std::vector<double>)> f,
    const double x,
    const std::vector<double> y,
    double h
); // rkf45

std::pair<std::vector<double>, std::vector<std::vector<double>>> driver
(
    const std::function<std::vector<double>(double, std::vector<double>)> f,
    const std::pair<double, double> interval,
    const std::vector<double> yi,
    double h=0.125,
    double acc=1e-3,
    double eps=1e-3
); // driver

std::vector<cubic_spline> interpolant
(
    const std::function<std::vector<double>(double, std::vector<double>)> f,
    const std::pair<double, double> interval,
    const std::vector<double> yi,
    double h=0.125,
    double acc=1e-3,
    double eps=1e-3
); // interpolant

#endif // RKF45_H

