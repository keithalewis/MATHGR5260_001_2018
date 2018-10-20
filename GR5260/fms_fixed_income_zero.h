// fms_fixed_income_zero.h - zero coupon bond
// C_u = 1
#pragma once
#include "fms_fixed_income_instrument.h"

namespace fms::fixed_income {

    template<class U = double, class C = double>
    class zero : public instrument<U,C> {
        U u;
        C c;
    public:
        zero(U u, C c = 1)
            : u(u), c(c)
        { }
    private:
        size_t _size() const override 
        {
            return 1;
        }
        const U* _time() const override 
        { 
            return &u; 
        }
        const C* _cash() const override 
        { 
            return &c; 
        }
    };

} // fms::fixed_income
