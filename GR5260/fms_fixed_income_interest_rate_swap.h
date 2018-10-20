// fms_fixed_income_interest_rate_swap.h - interest rates wap maturing at time u with coupon r and frequency f.
// t_j = j/f, j = 0, 1, ..., n = f*u
// C_j = r/f 0 < j < n
// C_0 = -1
// C_n = 1 + r*/f
#pragma once
#include <vector>
#include "fms_fixed_income_instrument.h"

namespace fms::fixed_income {

    template<class U = double, class C = double>
    class interest_rate_swap : public instrument<U,C> {
        std::vector<U> t;
        std::vector<C> c;
        frequency f;
    public:
        interest_rate_swap(U u, C r, frequency f)
            : t(/*???*/), c(/*???*/), f(f)
        {
            u = u;
            r = r;
            // t[j] = ???
            // c[j] = ???
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
