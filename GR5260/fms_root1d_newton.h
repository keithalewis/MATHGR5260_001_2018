// fms_root1d_newton.h - Newton's method 
#pragma once
#include <cmath>
#include <functional>
#include "fms_root1d.h"

namespace fms::root1d {

    template<class X = double, class Y = X, size_t max_iterations = 100>
    struct newton_solver : public abstract_solver<X> {
        X x;
        Y y;
        const std::function<Y(X)>& f;
        const std::function<Y(X)>& df;
        size_t n;
        newton_solver(X x, const std::function<Y(X)>& f, const std::function<Y(X)>& df)
            : x(x), f(f), df(df), n(0)
        { }
        newton_solver(const newton_solver&) = delete;
        newton_solver& operator=(const newton_solver&) = delete;
        virtual ~newton_solver()
        { }
        X _next()
        {
            ++n;
            y = f(x);
            x = x - y/df(x);

            return x;
        }
        bool _done() const 
        {
            if (n > max_iterations) {
                throw std::runtime_error("fms::root1d::newton_solver: exeeded maximum number of iterations");;
            }

            Y y_ = f(nextafter(x, X(1)));
            Y _y = f(nextafter(x, X(-1)));

            return fabs(y_) >= fabs(y) && fabs(_y) >= fabs(y);
        }
    };

} // fms::root1d