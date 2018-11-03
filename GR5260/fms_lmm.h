// fms_lmm.h - LIBOR Market Model
#pragma once
#include "fms_brownian.h"
/*
The LIBOR Market Model is parameterized by increasing times t_j,
futures quotes phi_j, at-the-money caplet volatilities, sigma_j, 
and a d x d correlation matrix, rho_{j,k}. 
The j-th future corresponds to the interval from t_{j-1} to t_j, j > 0
Just as for fsm::pwflat, we use the conventon t_{-1} = 0. We have
phi_0 is the cd rate and sigma_0 = 0.

Let Sigma Sigma' = [sigma_j rho_{j,k}] and
let Phi_j(t) = phi_j exp(Sigma B_t - ||Sigma||^2t/2) be the quote
at time t of the j-th futures, where B_t is d-dimensional standard
Brownian motion.

To account for convexity we let F_j(t) = Phi_j(t) - sigma_j^t (t_j - t)^2/2
be the forward at time t over the interval from t_j to t_{j+1}. In
our universal notation this is F_t(t_j, t_{j+1}).

Recall D(u) = exp(-int_0^u f(s) ds), where f(s) is the current forward curve, 
and D_t(u) = exp_t^u f_t(s) ds, where s -> f_t(s) is the forward curve at time t.
Given LMM data and a time t, we would like to generate a sample forward
curve at time t, s -> f_t(s). Note f_0(s) = f(s).


*/

namespace fms {

    template<class T = double, class F = double, class R = std::default_random_engine>
    class lmm {
        size_t n;
        std::vector<T> t;
        std::vector<F> phi;
        std::vector<F> sigma;
        std::brownian<F> B;
        R r;
    public:
        lmm(size_t n, const T* t, const F* phi, const F* sigma, const correlation<F>& e, R r))
            : n(n), t(t, t + n), phi(phi, phi + n), sigma(sigma, sigma + n), B(e), r(r)
        { }
        void reset()
        {
            B.reset();
        }
        // Populate t_ and f_ with sample forward curve at time u
        // return index of first t[i] > u
        size_t advance(T u, size_t n_, F* f_)
        {
            B.advance(u, r);

            //??? find first t[i] past u and generate random forwards using the LMM

            return 0; // ??? 
        }
    };

}