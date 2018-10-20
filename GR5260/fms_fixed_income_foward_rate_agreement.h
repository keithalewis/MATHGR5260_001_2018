// fms_fixed_income_forward_rate_agreement.h - forward rate agreement
// C_u = -1, C_v = 1 + f*(v - u)
#pragma once
#include <vector>
#include "fms_fixed_income_instrument.h"

namespace fms::fixed_income {

    template<class U = double, class C = double>
    class forward_rate_agreement : public instrument<U,C> {
        std::vector<U> t;
        std::vector<C> c;
    public:
        forward_rate_agreement(U u, U v, C f)
            : t(2), c(2)
        {
            t[0] = u;
            c[0] = -1;
            t[1] = v;
            c[1] = 1 + f*(u - v);
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
