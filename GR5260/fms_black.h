// fms_black.h - Black forward value pricing and greeks.
#pragma once
#include <cmath>

// The Black forward value of a put is E[max{k - F, 0}], where F is the forward and k is the strike.
// The forward is lognormal F = f exp(sigma B_t - sigma^2 t/2), where B_t is Brownian motion at time t.
//
// E[max{k - F, 0}] = E[(k - F) 1(F <= k)] = k P(F <= k) - E[F] P'(F <= k), 
// where dP'/dP = exp(s Z - s^2/2), s = sigma*sqrt(t), and Z is standard normal.
//
// Recall E[exp(r + sZ)] = exp(r + s^2/2), and E[exp(sZ) f(Z)] = E[f(Z + s)].
// The following are equivalent:
//   F <= k, sZ - s^2/2 <= log(k/f), Z <= s/2 + log(k/f)/s.
// P'(F <= k) = P'(Z <= s/2 + log(k/f)/s) = P(Z + s <= s/2 + log(k/f)/s) = P(Z <= -s/2 + log(k/f)/s)

namespace fms::black {

    // Standard normal cumulative distribution function.
    template<class X = double>
    inline X normal_cdf(X x) noexcept
    {
        static double sqrt2 = sqrt(X(2));

        return X(0.5) + erf(x / sqrt2) / 2;
    }

    template<class F = double, class S = double, class K = double>
    inline auto moneyness(F f, S s, K k)
    {
        return s/2 + log(k/f)/s;
    }

    // Black forward put value.
    template<class F = double, class S = double, class K = double>
    inline auto value(F f, S s, K k)
    {
        auto z = moneyness(f, s, k);

        return k * normal_cdf(z) - f * normal_cdf(z - s);
    }
    // Black forward put value with standard parameterization.
    template<class F = double, class S = double, class K = double, class T = double>
    inline auto value(F f, S sigma, K k, T t)
    {
        return value(f, sigma*sqrt(t), k);
    }

} // fms::black