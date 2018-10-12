// fms_black.h - Black forward value and greeks.
#pragma once
#include <cmath>
#include "fms_prob.h"
#include "fms_root1d_newton.h"
#include "../xll12/xll/ensure.h"

namespace fms::black {

    // The Black forward value of a put is E[max{k - F, 0}], 
    // where F is the forward and k is the strike.
    // The forward is lognormal F = f exp(-sigma^2 t/2 + sigma B_t), 
    // where B_t is Brownian motion at time t.
    //
    // E[max{k - F, 0}] = E[(k - F) 1(F <= k)] = k P(F <= k) - E[F] P'(F <= k), 
    // where dP'/dP = exp(s Z - s^2/2), s = sigma*sqrt(t), and Z is standard normal.
    //
    // Recall E[exp(r + sZ)] = exp(r + s^2/2), and E[exp(sZ - s^2/2) f(Z)] = E[f(Z + s)].
    // Note F <= k iff Z <= s/2 + log(k/f)/s = m, then moneyness.
    // P'(F <= k) = P'(Z <= m) = P(Z + s <= m) = P(Z <= m - s).

    // F <= k iff Z <= z
    template<class F = double, class S = double, class K = double>
    inline auto moneyness(F f, S s, K k)
    {
        return s/2 + log(k/f)/s;
    }

    // Black forward put value.
    template<class F = double, class S = double, class K = double>
    inline auto value(F f, S s, K k)
    {
        ensure(f >= 0);
        ensure(s >= 0);
        ensure(k >= 0);

        // If f = 0 then F = 0 so E max{k - F, 0} = k.
        // Note 1 + f == 1 is equivalent to fabs(f) < machine epsilon.
        if (1 + f == 1) {
            return k;
        }

        if (1 + k == 1) {
            return K(0);
        }

        if (1 + s == 0) {
            return std::max(k - f, F(0));
        }

        auto z = moneyness(f, s, k);

        return k * prob::normal_cdf(z) - f * prob::normal_cdf(z - s);
    }
    // Black forward put value with standard parameterization.
    template<class F = double, class S = double, class K = double, class T = double>
    inline auto value(F f, S sigma, K k, T t)
    {
        ensure(t >= 0);

        return value(f, sigma*sqrt(t), k);
    }

    // Derivative of value with respect to forward.
    // d/df E[max{k - F, 0}] = E[-exp(sZ - s^2/2)1(F <= k)] = -P'(F <= k)
    template<class F = double, class S = double, class K = double>
    inline auto delta(F f, S s, K k)
    {
        ensure(f >= 0);
        ensure(s >= 0);
        ensure(k >= 0);


        if (1 + k == 1) {
            return K(0);
        }
        // If f = 0 then F = 0 so E max{k - F, 0} = k.
        // Note 1 + f == 1 is equivalent to fabs(f) < machine epsilon.
        if (1 + f == 1) {
            return F(-1);
        }

        if (1 + s == 0) {
            return k == f ? S(0.5) : S(-1 * (k <= f));
        }

        auto z = moneyness(f, s, k);

        return  -prob::normal_cdf(z - s);
    }
    template<class F = double, class S = double, class K = double, class T = double>
    inline auto delta(F f, S sigma, K k, T t)
    {
        return delta(f, sigma*sqrt(t), k);
    }

    // Derivative of value with respect to volatility.
    template<class F, class S, class K, class T>
    inline auto vega(F f, S sigma, K k, T t)
    {
        ensure(f >= 0);
        ensure(sigma >= 0);
        ensure(k >= 0);
        ensure(t >= 0);

        // If f = 0 then F = 0 so E max{k - F, 0} = k.
        // Note 1 + f == 1 is equivalent to fabs(f) < machine epsilon.
        if (1 + f == 1) {
            return F(0);
        }

        // Test also for f == k???

        if (1 + k == 1) {
            return K(0);
        }

        if (1 + t == 0) {
            return T(0);
        }

        auto sqt = sqrt(t);
        auto s = sigma * sqt;
        auto z = moneyness(f, s, k);
        auto n = prob::normal_pdf(z);

        return f*n*sqt;
    }
    // Find Black put volatility with value v.
    template<class F, class V, class K, class T>
    inline auto implied(F f, V v, K k, T t)
    {
        V s0 = V(.2);
        std::function<V(V)> p = [f,v,k,t](V s) { return -v + value(f, s, k, t); }; 
        std::function<V(V)> dp = [f,k,t](V s) { return vega(f, s, k, t); };
        root1d::newton_solver<V, V> solver(s0, p, dp);

        return solver.solve(); 
    }

} // fms::black
