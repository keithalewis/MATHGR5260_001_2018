// fms_analytic.h - analytic numbers
// Analytic numbers allow one to compute derivatives to machine precision.
// Let J be an n x n matrix with J^{n-1} != 0 and J^n = 0.
// Analytic numbers have the form x = sum_{0 <= l < n} x_k J^k.
// For any analytic function f, f(xI + J) = f(x)I + f'(x)J + f''(x)/2 J^2 + ...
// This allows the derivatives of f to be computed using analytic numbers.
// We say n is the _order_ of the analytic number.
#pragma once
#include <algorithm>
#include <valarray>

namespace fms {
    
    // Toeplitz matrix where first row is valarray<X>
    template<class X>
    struct analytic {
        std::valarray<X> x;

        static constexpr X factorial(X n)
        {
            X n_ = 1; // n!

            while (n > 0) {
                n_ *= n;
                --n;
            };

            return n_;
        }
        static bool all_true(const std::valarray<bool>& x)
        {
            return std::all_of(std::begin(x), std::end(x), [](bool b) { return b;});
        }
        
        const std::valarray<X> left(size_t n) const
        {
            return x[std::slice(0, n, 1)];
        }
        std::valarray<X>& left(size_t n)
        {
            return x[std::slice(0, n, 1)];
        }
        
        const std::valarray<X> right(size_t n) const
        {
            return x[std::slice(order() - n, n, 1)];
        }
        std::valarray<X>& right(size_t n)
        {
            return x[std::slice(order() - n, n, 1)];
        }

        
        explicit analytic(size_t n)
            : x(n)
        { }
        // converting constructor for xI
        analytic(X x, size_t n)
            : x(n)
        {
            this->x[0] = x;
        }
        analytic(const analytic& y)
            : x(y.x)
        {
        }
        analytic& operator=(const analytic& y)
        {
            x.resize(y.x.size());
            x.operator=(y.x);

            return *this;
        }
        // yI
        analytic& operator=(X y)
        {
            x.operator=(X(0));
            x[0] = y;

            return *this;
        }
        ~analytic()
        { }

        bool operator==(const analytic& y) const
        {
            if (order() > y.order()) {
                return all_true(left(y.order()) == y.x)
                    && all_true(right(order() - y.order()) == X(0));
            }
            else if (order() < y.order()) {
                return all_true(x == y.left(order()))
                    && all_true(y.right(y.order() - order()) == X(0));
            }

            return all_true(x == y.x);
        }
        bool operator!=(const analytic& y) const
        {
            return !operator==(y);
        }

        size_t order() const
        {
            return x.size();
        }
        X operator[](size_t i) const
        {
            return x[i];
        }
        // i-th derivative
        X operator()(size_t i) const
        {
            return x[i]*factorial(i);
        }
        // results in x0 I
        analytic& resize(size_t n, const X x0 = X(0))
        {
            x.resize(n, X(0));
            x[0] = x0;

            return *this;
        }

        analytic& operator+=(const analytic& y)
        {
            if (order() > y.order()) {
                left(y.order()) += y.x;
            }
            else if (order() == y.order()) {
                x += y.x;
            }
            else {
                x += y.left(order());
            }

            return *this;
        }
        // non-member friend
        friend analytic operator+(analytic x, const analytic& y)
        {
            return x += y;
        }
    };
}
