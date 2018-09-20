// root1d.h - One dimensional root finding.
#pragma once
#include <functional>

namespace fms::root1d {

    // state, next, done
    template<class S, class N, class D>
    inline S solve(S s)
    {
        while (!D(s))
            s = N(s);

        return s;
    }

    template<class X = double>
    struct secant_state {
        X x, x_;
        std::function<X(X)> F;
        secant_state(X x, X x_, const std::function<X(X)>& F)
            : x(x), x_(x_), F(F)
        { }
        secant_state& next()
        {
            X f = F(x);
            X f_ = F(x_);

            x = (x*f_ - x_*f)/(f_ - f);
            std::swap(x, x_);

            return *this;
        }
    };
    template<class X = double>
    inline secant_state<X> secant_next(secant_state<X> s)
    {
        return s.next();
    }

    template<class X = double>
    struct secant_done {
        secant_done()
        { }
    };


} // fms::root1d
