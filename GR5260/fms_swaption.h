// fms_swaption.h - Swaption pricing
#pragma once
#include <functional>
#include "fms_fixed_income_interest_rate_swap.h"
#include "fms_lmm.h"

namespace fms {

    // Option on a swap paying
    // C_t = max{F_t(t_0,...,t_n) - k, 0}
    template<class T = double, class F = double>
    F swaption_payoff(T tenor, fixed_income::frequency freq, F k, T t, const std::function<F(T)>& D)
    {
        
        auto Dt = D(t);
        // Discount curve at time t.
        std::function<F(T)> Du = [Dt,&D](T u) { return D(u)/Dt; };
        F c = fixed_income::interest_rate_swap<T,F>::par_coupon(tenor, freq, Du);
     
        return std::max(c - k, F(0));
    }
    template<class T = double, class F = double>
    F swaption(T tenor, fixed_income::frequency freq, F k, T t, lmm<T,F>& lmm, size_t N = 10000)
    {
        std::default_random_engine dre; // should be passed in
        F pv = F(0);

        for (size_t n = 1; n <= N; ++n) {
            std::vector<F> f(lmm.size()); // foward curve
            
            // Initial forward curve
            lmm.reset();
            lmm.advance(0, f.data(), dre);

            // Deterministic discount to first futures settlement.
            F Dt = pwflat::discount(lmm.t[0], lmm.size(), lmm.t.data(), f.data());

            // Generate (approximate) stochastic discount to time t
            for (size_t j = 0; j < lmm.size() && lmm.t[j] < t; ++j) {
                lmm.advance(lmm.t[j], f.data(), dre);
                Dt *= exp(-f[j]*(lmm.t[j+1] - lmm.t[j]));
            }
            
            size_t k = lmm.advance(t, f.data(), dre);
            Dt *= exp(-f[k]*(t - lmm.t[k]));

            std::function<F(T)> D = [&f, &lmm](T t) { 
                return pwflat::discount(t, lmm.size(), lmm.t.data(), f.data());
            };
            
            pv += (swaption_payoff<T,F>(tenor, freq, k, t, D) - pv)/n;
        }

        return pv;
    }

}