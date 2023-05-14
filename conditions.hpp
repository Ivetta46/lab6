#pragma once

#include <cmath>

// Here defind functions for initial and boundary conditions.

constexpr double tau = 0.2;
constexpr double hx = 0.2;
constexpr double hy = 0.2;
constexpr double epsilon = 0.001;

double left(double y, double t)
{
    return std::pow(1 + y, 4.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double right(double y, double t)
{
    return std::pow(2 + y, 4.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double bottom(double x, double t)
{
    return std::pow(1 + x, 4.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double top(double x, double t)
{
    return std::pow(2 + x, 4.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double initial(double x, double y)
{
    return std::pow(1 + x + y, 4.0 / 3.0) / std::pow(100, 1.0 / 3.0);
}

double K(double arg)
{
    return std::pow(arg, 1.5);
}