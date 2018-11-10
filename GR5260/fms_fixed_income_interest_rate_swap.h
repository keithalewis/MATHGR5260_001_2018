// fms_fixed_income_interest_rate_swap.h - interest rates wap maturing at time u with coupon r and frequency f.
// t_j = j/f, j = 0, 1, ..., n = f*u
// C_0 = -1
// C_j = r/f, 0 < j < n
// C_n = 1 + r*/f
#pragma once
#include <vector>
#include "fms_fixed_income_instrument.h"

namespace fms::fixed_income {

    template<class U = double, class C = double>
    class interest_rate_swap : public instrument<U,C> {
        U u;
        C r;
        frequency f;
        std::vector<U> t;
        std::vector<C> c;
    public:
        interest_rate_swap(U u, C r, frequency f)
            : u(u), r(r), f(f), t(1 + static_cast<U>(f)*u), c(1 + static_cast<U>(f)*u)
        {
            t[0] = 0;
            c[0] = -1;
            for (size_t i = 1; i < t.size(); ++i) {
                t[i] = i/static_cast<U>(f);
                c[i] = r/static_cast<U>(f);
            }
            c.back() += 1;
        }
        // F^delta(t_0,...,t_n) = (1 - D(t_n)/sum_1^n delta_j D(t_j)
        // ??? Assume the daycount fraction is delta_j = t[j-1] - t[j].
        C par_coupon(const std::function<C(U)>& D)
        {
            C c = C(0);

            //??? Calculate the par coupon assuming D is the discount.

            return c;
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

} // fms::fixed_income
