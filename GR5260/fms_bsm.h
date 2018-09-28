// fms_bsm.h - Black-Scholes/Merton valuation and greeks.
#pragma once
#include "fms_black.h"

namespace fms {
    namespace bsm {

        // The Black-Scholes/Merton value of a put is exp(-rt)E[max{k - S, 0}], 
        // where S is the stock price and k is the strike.
        // The stock is lognormal S = s exp((r - sigma^2/2)t + sigma B_t), 
        // where B_t is Brownian motion at time t.
        // The relationship between the B-S/M value and the Black value is
        // v_BMS(r, s, sigma, k, t) = exp(-rt) v_B(s e^rt, sigma, k, t); 
        //
        // Black-Scholes/Merton put value.
        template<class R = double, class F = double, class S = double, class K = double, class T = double>
        inline auto value(R r, F s, S sigma, K k, T t)
        {
            auto D = exp(-r * t);

            return D * black::value(s / D, sigma, k, t);
        }

        // d/ds v_BMS(r, s, sigma, k, t) 
        //  = d/ds exp(-rt) v_B(s exp(rt), sigma, k, t)
        //  = exp(-rt) d/df v_B(s exp(rt), sigma, k, t) exp(rt)
        //  = d/df(s exp(rt), sigma, k, t)
        // Black-Scholes/Merton put delta.
        template<class R = double, class F = double, class S = double, class K = double, class T = double>
        inline auto delta(R r, F s, S sigma, K k, T t)
        {
            auto f = s * exp(r*t);

            return black::delta(f, sigma, k, t);
        }

    } // bsm
} // fms