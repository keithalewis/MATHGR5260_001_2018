// fms_analytic.h - analytic numbers
#pragma once
#include <valarray>

namespace fms {
    
    // Toeplitz matrix where first row is valarray<X>
    template<class X>
    struct analytic {
        std::valarray<X> x;
        explicit analytic(size_t n)
            : x(n)
        { }
        // converting constructor for xI
        analytic(X x, size_t n)
            : x(n)
        {
            this->x[0] = x;
        }
        analytic(const analytic& x)
            : x(x)
        {
        }
        analytic operator=(const analytic& x)
        {
            this->x.resize(x.size());
            this->x.operator=(x);
        }
        ~analytic()
        { }

        size_t size() const
        {
            return x.size();
        }
        analytic& resize(size_t n, const X x0 = X(0))
        {
            x.resize(n, x0);

            return *this;
        }

        analytic& operator+=(const analytic& y)
        {
            x += y.x;

            return *this;
        }
        // non-member friend
        friend analytic operator+(analytic x, const analytic& y)
        {
            return x += y;
        }
    };
}
