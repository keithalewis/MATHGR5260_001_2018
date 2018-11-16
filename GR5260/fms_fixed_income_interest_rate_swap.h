// fms_fixed_income_interest_rate_swap.h - interest rates wap maturing at time u with coupon r and frequency f.
// t_j = j/f, j = 0, 1, ..., n = f*u
// C_0 = -1
// C_j = r/f, 0 < j < n
// C_n = 1 + r*/f
#pragma once
#include <functional>
#include <vector>
#include "fms_fixed_income_instrument.h"

namespace fms::fixed_income {

    template<class U = double, class C = double>
    class interest_rate_swap : public instrument<U,C> {
        U u;
        C r;
        frequency q;
        std::vector<U> t;
        std::vector<C> c;
    public:
        interest_rate_swap(U u, C r, frequency q)
            : u(u), r(r), q(q), t(static_cast<size_t>(1 + static_cast<U>(q)*u)), c(static_cast<size_t>(1 + static_cast<U>(q)*u))
        {
            U dt = static_cast<U>(q);
            dt = 1/dt;
            t[0] = 0;
            c[0] = -1;
            for (size_t i = 1; i < t.size(); ++i) {
                t[i] = i*dt;
                c[i] = r*dt;
            }
            c.back() += 1;
        }
    private:
        size_t _size() const override 
        {
            return t.size();
        }
        const U* _time() const override 
        { 
            return t.data(); 
        }
        const C* _cash() const override 
        { 
            return c.data(); 
        }
    };

    // F^delta(t_0,...,t_n) = (D(t_0) - D(t_n))/sum_1^n delta_j D(t_j)
    template<class U = double, class C = double>
    inline C par_coupon(const interest_rate_swap<U,C>& irs, const std::function<C(U)>& D)
    {
        C sum = C(0);
        const U* u = irs.time();

        for (size_t j = 1; j < irs.size(); ++j) {
            sum += (u[j] - u[j - 1])*D(u[j]);
        }

        return (D(u[0]) - D(u[irs.size()- 1]))/sum;
    }

} // fms::fixed_income
