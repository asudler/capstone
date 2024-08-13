#ifndef RKF45_H
#define RKF45_H

#include <functional>
#include <utility>
#include <vector>
#include "/home/asudler/git/capstone/spline/spline.h"

template <typename T>
std::pair<std::vector<T>, std::vector<T>> rkf45
(
    const std::function<std::vector<T>(double, std::vector<T>)> f,
    const double x,
    const std::vector<T> y,
    double h
); // rkf45

template <typename T>
std::pair<std::vector<double>, std::vector<std::vector<T>>> driver
(
    const std::function<std::vector<T>(double, std::vector<T>)> f,
    const std::pair<double, double> interval,
    const std::vector<T> yi,
    double h=0.125,
    const double acc=1e-3,
    const double eps=1e-3
); // driver

template <typename T>
std::vector<cubic_spline<T>> interpolant
(
    const std::function<std::vector<T>(double, std::vector<T>)> f,
    const std::pair<double, double> interval,
    const std::vector<T> yi,
    double h=0.125,
    const double acc=1e-3,
    const double eps=1e-3
); // interpolant

#endif // RKF45_H

