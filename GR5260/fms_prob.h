// fms_prob.h - Probability related functions
#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

namespace fms::prob {

    // Standard normal probability density function.
    template<class X = double>
    inline X normal_pdf(X x) noexcept
    {
        static X sqrt2pi = sqrt(X(2)*X(M_PI));

        return exp(-x*x/2)/sqrt2pi;
    }
    // Standard normal cumulative distribution function.
    template<class X = double>
    inline X normal_cdf(X x) noexcept
    {
        static X sqrt2 = sqrt(X(2));

        return X(0.5) + erf(x/sqrt2)/2;
    }

} // fms::prob