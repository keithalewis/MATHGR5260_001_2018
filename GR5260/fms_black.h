// fms_black.h - Black forward value pricing and greeks.
#pragma once
#include <cmath>

namespace fms {

    // Standard normal cumulative distribution function.
    template<class X = double>
    inline X normal_cdf(X x) noexcept
    {
        static X sqrt2 = sqrt(X(2));

        return X(0.5) + erf(x / sqrt2) / 2;
    }

    namespace black {

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

        // d/dF E[max{k - F, 0}] = E[-exp(sZ - s^2/2)1(F <= k)] = -P'(F <= k)
        // Black forward put delta
        template<class F = double, class S = double, class K = double>
        inline auto delta(F f, S s, K k)
        {
            return f*s*k; // ???Implement the formula above
        }
        template<class F = double, class S = double, class K = double, class T = double>
        inline auto delta(F f, S sigma, K k, T t)
        {
            return delta(f, sigma*sqrt(t), k);
        }
    } // black
    namespace bsm {

        // The Black-Scholes/Merton value of a put is exp(-rt)E[max{k - S, 0}], 
        // where S is the stock price and k is the strike.
        // The stock is lognormal S = s exp((r - sigma^2/2)t + sigma B_t), 
        // where B_t is Brownian motion at time t.
        // The relationship between the B-S/M value and the Black value is
        // v_BMS(r, s, sigma, k, t) = exp(-rt) v_B(s e^rt, sigma, k, t); 
        //
        // Black-Scholes/Merton put value.
        template<class F = double, class S = double, class K = double, class T = double>
        inline auto value(F f, S sigma, K k, T t)
        {
            return 0; //???Implement in terms of the black::value
        }

        // d/ds v_BMS(r, s, sigma, k, t) = ???Insert formula here
        // Black-Scholes/Merton put delta.
        template<class F = double, class S = double, class K = double, class T = double>
        inline auto delta(F f, S sigma, K k, T t)
        {
            return 0; //???Implement in terms of the black::delta
        }

    } // bsm
} // fms