// fms_ho_lee.h - Ho-Lee normal short rate model
#pragma once
#include "fms_bsm.h"

/* ???
The Ho-Lee model for the short rate is

    f_t = phi(t) + sigma B_t,

where B_t is standard Brownian motion.

The _stochastic discount_ is D_t = exp(-int_0^t f_s ds) = exp(-int_0^t [phi(s) + sigma B_s] ds)

The price at time t of a zero coupon bond maturing at time u is

D_t(u) = E_t D_u/D_s = E_t exp(-int_t^u phi(s) ds + int_t^u sigma B_s ds)
                     = exp(-int_t^u phi(s) ds) E_t exp(-int_t^u sigma B_s ds)

Since d(s B_s) = s dB_s + B_s ds (+ ds dB_s = 0), we have

    int_t^u B_s ds = (u B_u - t B_t) - int_t^u s dB_s
                   = (u B_u - u B_t + u B_t - t B_t) - int_t^u s dB_s
                   = (u - t)B_t + int_t^u (u - s) dB_s.

So

    E_t exp(-int_t^u sigma B_s ds) = E_t exp(-sigma(u - t)B_t + int_t^u sigma(u - s) dB_s)
                                   = exp(-sigma(u - t)B_t + int_t^u sigma^2(u - s)^2 ds/2),

using the fact that exp(int_0^t a(s) dB_s - int_0^t a(s)^2 ds/2) is a martingale.

This shows

    D_t(u) = exp(-sigma(u - t)B_t - int_t^u [phi(s) - sigma^2(u - s)^2/2] ds).

In particular

    D(t) = D_0(t)
         = exp(-int_0^t [phi(s) - sigma^2(t - s)^2/2] ds)
         = exp(-int_0^t phi(s) ds + int_0^t sigma^2 (t - s)^2/2 ds)
         = exp(-int_0^t phi(s) ds + int_0^t sigma^2 s^2/2 ds)
         = exp(-int_0^t [phi(s) - sigma^2 s^2/2] ds)

Since D(t) = exp(-int_0^t f(s) ds) we have f(t) = phi(t) - sigma^2 t^2/2.

A (interest rate) _floor_ over the interval [u,v] with strike k pays
max{k - F^delta_u(u,v), 0} at time v, where F^delta_u(u,v) is the forward rate
at time u over the interval [u,v] using day count basis delta. 

Recall the forward at time t over [u,v] is F^delta_t(u,v) = (D_t(u)/D_t(v) - 1)/dcf
where dcf is the day count fraction.

The floor value is

    E max{k - F,0} D_v = E max{k - (1/D_u(v) - 1)/delta, 0} D_v
                       = ???
                       = ??? Document the derivation of the formula below here.
                       = ???
                       = D(u) E max{(k + 1/delta)D_u(v)e^gamma - 1/delta, 0}

where gamma = Cov(log D_u(v), log D_u)
                       
*/

namespace fms::ho_lee {

    // Expected value of log D_t(u) = -int_t^u [phi(s) - sigma^2(u - s)^2/2] ds
    template<class X = double>
    inline auto ElogD_(X t, X u, X Dt, X Du, X sigma)
    {
        return 0; //???
    }

    // Variance of log D_t(u) = -int_t^u [phi(s) - sigma^2(u - s)^2/2] ds
    template<class X = double>
    inline auto VarlogD_(X t, X u, X Dt, X Du, X sigma)
    {
        return 0; //???
    }

    // Covariance of log D_u(v) and log D_u.
    template<class X = double>
    inline auto CovlogD_(X u, X v, ...)
    {
        return 0; //???
    }

    // E max{k - F,0} D_v = D(u) E max{(k + 1/delta)D_u(v)e^gamma - 1/delta, 0}
    // Use Black-Scholes/Merton to value an interest rate floor.
    template<class X = double>
    inline auto floor(X k, X delta, X u, X v, X Du, X Dv, X sigma)
    {
        X R = -log(Du)/t; 
        X S = 0; //??? B-S/M spot
        X Sigma = 0; //??? B-S/M volatility
        X K = 0; //??? B-S/M strike

        // Black-Scholes/Merton call value.
        return S - K*Du + bms::value(R, S, Sigma, K, u);
    }
}