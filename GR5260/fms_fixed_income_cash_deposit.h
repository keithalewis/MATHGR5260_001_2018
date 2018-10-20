// fms_fixed_income_cash_deposit.h - cash deposit
// C_u = 1 + r*u
#pragma once
#include "fms_fixed_income.h"

namespace fms::fixed_income {

    template<class U = double, class C = double>
    class cash_deposit : public instrument<U,C> {
        U u;
        C c;
    public:
        cash_deposit(U u, C r = 0)
            : u(u), c(1 + r*u)
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
