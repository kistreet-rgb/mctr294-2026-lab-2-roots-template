#include "roots.hpp"
#include <cmath>
#include <algorithm> // For std::abs

// Constants defined by lab requirements [cite: 26, 27]
static const double TOLERANCE = 1e-6;
static const int MAX_ITER = 1000000;

bool bisection(std::function<double(double)> f, double a, double b, double *root)
{
    // Check if interval is valid (signs must be opposite)
    if (f(a) * f(b) >= 0) // the inpuut a and b must have diffent signs otherwise the algorithm will not work.

    {
        return false;
    }

    double c = a;
    for (int i = 0; i < MAX_ITER; ++i)
    {
        // Calculate midpoint
        c = (a + b) / 2.0;
        double fc = f(c);

        // Check for convergence
        // we check if f(c) is close to zero or if the interval is smaller than the tolderance that we cant seach for the root any longer

        if (std::abs(fc) < TOLERANCE || std::abs(b - a) < TOLERANCE)
        {
            *root = c;
            return true;
        }

        // Update bracket
        if (fc * f(a) < 0)// weneed to check if the root is in the left or right side of the input.
        {
            b = c; // Root is in left half
        }
        else
        {
            a = c; // Root is in right half
        }
    }

    return false; // Failed to converge within max iterations we definned earlier 
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root)//uses the xint of the line to approximate the root its both simplee and effecient 
{
    // Check if interval is valid (signs must be opposite)
    if (f(a) * f(b) >= 0)
    {
        return false;
    }

    double c = a;
    for (int i = 0; i < MAX_ITER; ++i)
    {
        double fa = f(a);
        double fb = f(b);

        // Prevent division by zero
        if (std::abs(fb - fa) < 1e-15)
        {
            return false;
        }

        // False Position formula
        c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        // Check for convergence
        if (std::abs(fc) < TOLERANCE)
        {
            *root = c;
            return true;
        }

        // Update bracket similar to the bisection method we have to check what side the root is on.
        if (fc * fa < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }

    return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) //we say g is a derivative of f we use the derivative to find the root faster but it requries a good initial guess 
{
    double x = c;

    // Check if initial guess is valid
    if (x < a || x > b)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double fx = f(x);
        double gx = g(x);

        // Check success
        if (std::abs(fx) < TOLERANCE)
        {
            *root = x;
            return true;
        }

        // Check derivative to avoid divide by zero
        if (std::abs(gx) < 1e-15) // we want the derivative to not be zero because we cant divide by zero also if the derivative is too small the method wont converge
        {
            return false;
        }

        // Newton step: x_new = x - f(x)/f'(x)
        double x_new = x - (fx / gx);

        // Check bounds 
        // "Returns false if ... iteration leaves the interval"
        if (x_new < a || x_new > b)
        {
            return false;
        }

        x = x_new;
    }

    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    // Secant requires two points. We use c (guess) as x0.
    // We generate x1 by adding a small perturbation to c.
    double x0 = c;
    double x1 = c + 0.001;

    // Ensure x1 is within bounds. If c is too close to b, go backwards.
    if (x1 > b)
    {
        x1 = c - 0.001;
    }
    
    // Initial bounds check
    if (x0 < a || x0 > b || x1 < a || x1 > b) return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double f0 = f(x0);
        double f1 = f(x1);

        if (std::abs(f1) < TOLERANCE)
        {
            *root = x1;
            return true;
        }

        // Avoid divide by zero
        if (std::abs(f1 - f0) < 1e-15)
        {
            return false;
        }

        // Secant formula: x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        // Check bounds to see if the new point is within the given interval
        if (x2 < a || x2 > b)
        {
            return false;
        }

        // Shift th points for next iteration 
        x0 = x1;
        x1 = x2;
    }

    return false;
}
//end of code! it was created by Keiron Street with refrences from  AI for structures and W3Schools.com for syntax help.
// thank you for using my code!
