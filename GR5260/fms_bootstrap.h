// fms_bootstrap.h - Bootstrap a piecewise flat forward curve.
#pragma once
#include <limits>
#include "fms_root1d_newton.h"
#include "fms_pwflat.h"

namespace fms::pwflat {
	
	// Extrapolate curve to match price with present value.
	template<class T, class F>
	inline std::pair<T,F> bootstrap(F p, 
        size_t m, const T* u, const F* c, 
        size_t n, const T* t, const F* f, F _f = 0)
	{
        // expiration must be past the end of the forward curve
        ensure (m > 0);
        ensure (n == 0 || u[m-1] > t[n-1]);
        
        // end of curve 
        auto t_ = n == 0 ? 0 : t[n-1];
        // discount to end
        auto D_ = discount(t_, n, t, f);
        // last cash flow
        auto c_ = c[m - 1];
        // last cash flow time
        auto u_ = u[m - 1];

        // If only one cash flow occurs past the end of the curve there is a closed form solution:
        // We have p = pv + c D e^{-f(u - t)}, where pv is the present value of all but the last
        // cash flow, c is the last cash flow, and u is the last cash flow time.
        if (m == 1 || u[m - 2] <= t_) {
            auto pv = present_value(m - 1, u, c, n, t, f);
            
            return std::make_pair(u_, log((p - pv)/(c_*D_))/(t_ - u_));
        }

        // If exactly two cash flows and price is 0, then we know u[0] > t_ or else
        // the previous case would hold.
        // 0 = c0 D exp(-f(u0 - t)) + c1 D exp(-f(u1 -t)) so
        // f = - log(-c0/c1)/(u1 - u0).
        if (p == 0 && m == 2) {
            ensure (u[0] > t_);
            ensure (u[0] < u[1]);

            return std::make_pair(u_, log(-c[0]/c[1])/(u[0] - u[1]));
        }

        std::function<F(F)> pv = [p, m, u, c, n, t, f](F _f) {
			return -p + present_value(m, u, c, n, t, f, _f);
		};
		std::function<F(F)> dpv = [m, u, c, n, t, f](F _f) {
			return partial_duration(m, u, c, n, t, f, _f);
		};

		// Use fms::root1d::newton to solve for the extrapolated value.
		if (n > 0)
			_f = f[n-1];

        //???Check if a solution is possible.

        root1d::newton_solver<F,F> solver(pv, dpv);

		_f = solver.solve(_f);

		return std::make_pair(u_,_f);
    }
}
