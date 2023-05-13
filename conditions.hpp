#pragma once

#include <cmath>

constexpr double tau = 1;
constexpr double hx = 1;
constexpr double hy = 1;
constexpr double epsilon = 0.1;

void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
    {
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}

double bottom(double y, double t)
{
    return std::pow(1 + y, 5.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double top(double y, double t)
{
    return std::pow(2 + y, 5.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double left(double x, double t)
{
    return std::pow(1 + x, 5.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double right(double x, double t)
{
    return std::pow(2 + x, 5.0 / 3.0) / std::pow(10 - 28 * t / 3.0, 2.0 / 3.0);
}

double initial(double x, double y)
{
    return std::pow(1 + x + y, 4.0 / 3.0) / std::pow(100, 1.0 / 3.0);
}

double K(double arg)
{
    return std::pow(arg, 1.5);
}