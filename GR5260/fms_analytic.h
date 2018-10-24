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
    class analytic {
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
        auto left(size_t n) const // take
        {
            return std::slice(0, n, 1);
        }        
        auto right(size_t n) const
        {
            return std::slice(order() - n, n, 1);
        }
    public:
        explicit analytic(size_t n)
            : x(n)
        { }
        // converting constructor for x0 I + J
        analytic(X x0, size_t n)
            : x(n)
        {
            if (n > 0)
                x[0] = x0;
            if (n > 1)
                x[1] = X(1);
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
                return (x[left(y.order())] == y.x).min() == true
                    && (x[right(y.order())] == X(0)).min() == true;
            }
            else if (order() < y.order()) {
                return (x == y.x[left(order())]).min() == true
                    && (y.x[right(y.order())] == X(0)).min() == true;
            }

            return (x == y.x).min() == true;
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
                x[left(y.order())] += y.x;
            }
            else if (order() < y.order()) {
                x += y.x[left(order())];
            }
            else {
                x += y.x;
            }

            return *this;
        }
        // non-member friend
        friend analytic operator+(analytic x, const analytic& y)
        {
            return x += y;
        }
        analytic& operator-=(const analytic& y)
        {
            if (order() > y.order()) {
                x[left(y.order())] -= y.x;
            }
            else if (order() < y.order()) {
                x -= y.x[left(order())];
            }
            else {
                x -= y.x;
            }

            return *this;
        }
        // non-member friend
        friend analytic operator-(analytic x, const analytic& y)
        {
            return x -= y;
        }
        analytic& operator*=(const analytic& y)
        {
            std::valarray<X> z(order());
            for (size_t i = 0; i < z.size(); ++i) {
                for (size_t j = 0; j < x.size() && x.size() - j < y.order(); ++j) {
                    z[i] += x[j]*y.x[z.size() - j];
                }
            }
            std::swap(x, z);

            return *this;
        }
        // non-member friend
        friend analytic operator*(analytic x, const analytic& y)
        {
            return x *= y;
        }

    };
}
