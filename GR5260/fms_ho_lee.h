// fms_ho_lee.h - Ho-Lee normal short rate model
#pragma once
#include "fms_bsm.h"

/*
The Ho-Lee model for the short rate is

    f_t = phi(t) + sigma B_t,

where B_t is standard Brownian motion.

The _stochastic discount_ is D_t = exp(-int_0^t f_s ds) = exp(-int_0^t [phi(s) + sigma B_s] ds)
so int_0^t f_s ds = int_0^t [phi(s) + sigma B_s] ds
and E int_0^t f_s ds = int_0^t phi(s) ds
hence E f_t = phi(t).

The price at time t of a zero coupon bond maturing at time u is

D_t(u) = E_t D_u/D_s = E_t exp(-int_t^u phi(s) ds + int_t^u sigma B_s ds)
                     = exp(-int_t^u phi(s) ds) E_t exp(-int_t^u sigma B_s ds)

Since d(s B_s) = s dB_s + B_s ds (+ ds dB_s = 0), we have

    int_t^u B_s ds = (u B_u - t B_t) - int_t^u s dB_s
                   = (u B_u - u B_t + u B_t - t B_t) - int_t^u s dB_s
                   = (u - t)B_t + int_t^u (u - s) dB_s.

So

    E_t exp(-int_t^u sigma B_s ds) = E_t exp(-sigma(u - t)B_t - int_t^u sigma(u - s) dB_s)
                                   = exp(-sigma(u - t)B_t + int_t^u sigma^2(u - s)^2/2 ds),

because exp(int_0^t a(s) dB_s - int_0^t a(s)^2 ds/2) is a martingale
and so E_t^u exp(int_t^ a(s) dB_s) = exp(int_t^u a(s)^2/2 ds).

This shows

    D_t(u) = exp(-sigma(u - t)B_t - int_t^u [phi(s) - sigma^2(u - s)^2/2] ds).

In particular

    D(t) = D_0(t)
         = exp(-int_0^t [phi(s) - sigma^2(t - s)^2/2] ds)
         = exp(-int_0^t phi(s) ds + int_0^t sigma^2 (t - s)^2/2 ds)
         = exp(-int_0^t phi(s) ds + int_0^t sigma^2 s^2/2 ds)
         = exp(-int_0^t [phi(s) - sigma^2 s^2/2] ds)

Since D(t) = exp(-int_0^t f(s) ds) we have f(t) = phi(t) - sigma^2 t^2/2.

Using E exp(N) = exp(E[N] + Var[N]/2) and 

    Var(int_0^t B_s ds) = Cov(int_0^t B_u du, int_0^t B_v dv) 
                        = int_0^t int_0^t Cov(u,v) du dv
                        = int_0^t int_0^t min{u,v} du dv
                        = 2 int_0^t int_0^v u du dv
                        = 2 int_0^t v^2/2 dv
                        = t^3/3

we have D(t) = exp(-int_0^t phi(s) ds + sigma^2 t^3/6). This is the same as
the formula above since int_0^t s^2/2 ds = t^3/6.

Note int_t^u (u - s)^2/2 ds = -(u - s)^3/6|_t^u = (u - t)^3/6
and  int_t^u s^2/2 ds = (u^3 - t^3)/6.
Since (u^3 - t^3) - (u - t)^3 = 3 u^2 t - 3 u t^2 we get
    
    D_t(u) = exp(-sigma(u - t)B_t) D(u)/D(t) exp(-sigma^2 ut(u - t)/2)
*/

namespace fms::ho_lee {

    // Expected value of log D_t(u) = log D(u)/D(t) - sigma^2 ut(u - t)/2.
    template<class X = double>
    inline auto ElogD_(X t, X u, X Dt, X Du, X sigma)
    {
        return log(Du/Dt) - sigma*sigma*u*t*(u - t)/2;
    }

    // Variance of log D_t(u) = sigma^2(u - t)^2 t
    template<class X = double>
    inline auto VarlogD_(X t, X u, X sigma)
    {
        return sigma*sigma*(u - t)*(u - t)*t;
    }

    // The covariance of log D_t(u) and log D_t is
    // Cov(sigma(u - t)B_t, int_0^u sigma B_s ds)
    // = sigma^2 (u - t) [int_0^u Cov(B_t,B_s) ds]
    // = sigma^2 (u - t) [int_0^u min{t,s} ds]
    // = sigma^2 (u - t) [int_0^t s ds + int_t^u t ds]
    // = sigma^2 (u - t) [t^2/2 + t(u - t)]
    // = sigma^2 (u - t) t (t/2 + (u - t))
    template<class X = double>
    inline auto CovlogD_(X t, X u, X sigma)
    {
        return sigma*sigma*(u - t)*t*(u - t/2);
    }

    /*
    A (interest rate) _floor_ over the interval [u,v] with strike k pays
    max{k - F^delta_u(u,v), 0} at time v, where F^delta_u(u,v) is the forward rate
    at time u over the interval [u,v] using day count basis delta. 

    Recall the forward at time t over [u,v] is F^delta_t(u,v) = (D_t(u)/D_t(v) - 1)/dcf
    where dcf is the day count fraction corresponding to delta, u, and v;

    The floorlet value is

        E max{k - F,0} D_v = E max{k - (1/D_u(v) - 1)/dcf, 0} D_v
                           = E E_u max{k - (1/D_u(v) - 1)/dcf, 0} D_v
                           = E max{k - (1/D_u(v) - 1)/dcf, 0} E_u D_v
                           = E max{k - (1/D_u(v) - 1)/dcf, 0} D_u(v) D_u
                           = E max{kD_u(v) - (1 - D_u(v))/dcf, 0} D_u
                           = E max{(k + 1/dcf)D_u(v) - 1/dcf, 0} D_u
                           = D(u) E max{(k + 1/dcf)D_u(v)e^gamma - 1/dcf, 0}
        
    where gamma = Cov(log D_u(v), log D_u).
    */

    // E max{k - F,0} D_v = D(u) E max{(k + 1/dcf)D_u(v)e^gamma - 1/dcf, 0}
    // Use Black-Scholes/Merton to value an interest rate floor.
    // D(u) = exp(-R u), S_u = (k + 1/dcf)D_u(v)e^gamma and k = 1/dcf.
    template<class X = double>
    inline auto floor(X k, X dcf, X u, X v, X Du, X Dv, X sigma)
    {
        X R = -log(Du)/u;
        X Elog = ElogD_(u, v, Du, Dv, sigma);
        X Varlog = VarlogD_(u, v, sigma);
        X gamma = CovlogD_(u, v, sigma);
        X S = (k + 1/dcf)*exp(Elog + Varlog/2 + gamma)/Du;
        X Sigma = sqrt(Varlog/u); // Sigma^2 u = Varlog
        X K = 1/dcf;

        // Black-Scholes/Merton call value.
        return S - K*Du + bsm::value(R, S, Sigma, K, u);
    }
}