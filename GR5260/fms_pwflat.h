// fms_pwflat.h - piecewise flat curve
/*
    f(t) = f[i] if t[i-1] < t <= t[i]
         = _f   if t > t[n-1]
    and undefined if t < 0
    Note f(0) = f[0].

    |                                   _f
    |        f[1]             f[n-1] (--------
    | f[0] (----- ...       (------]
    [------]      ... ------]
    |
    0-----t[0]--- ... ---t[n-2]---t[n-1]
*/
#pragma once
#include <cmath>     // exp
#include <algorithm> // adjacent_find
#include <limits>    // quiet_Nan()
#include <numeric>   // upper/lower_bound

namespace fms::pwflat {

    // strictly increasing values
    template<class T>
    inline bool strictly_increasing(size_t n, const T* t) noexcept
    {
        return t + n == std::adjacent_find(t, t + n, std::greater_equal<T>{});
    }

    // piecewise flat curve
    // return f[i] if t[i-1] < u <= t[i], _f if u > t[n-1]
    // assumes t[i] monotonically increasing
    template<class T, class F>
    inline F value(const T& u, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        if (u < 0 || !strictly_increasing(n, t))
            return std::numeric_limits<F>::quiet_NaN();

        if (n == 0)
            return _f;

        auto ti = std::lower_bound(t, t + n, u);

        return ti == t + n ? _f : f[ti - t];
    }

    // int_0^u f(t) dt
    template<class T, class F>
    inline F integral(const T& u, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        if (u < 0 || !strictly_increasing(n, t))
            return std::numeric_limits<F>::quiet_NaN();

        F I{ 0 };
        T t_{ 0 };

        size_t i;
        for (i = 0; i < n && t[i] <= u; ++i) {
            I += f[i] * (t[i] - t_);
            t_ = t[i];
        }
        I += (n == 0 || u > t[n - 1] ? _f : f[i]) *(u - t_);

        return I;
    }

    // discount D(u) = exp(-int_0^u f(t) dt)
    template<class T, class F>
    inline F discount(const T& u, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        return exp(-integral(u, n, t, f, _f));
    }

    // spot r(u) = (int_0^u f(t) dt)/u
    template<class T, class F>
    inline F spot(const T& u, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        return u <= t[0] ? f[0] : integral(u, n, t, f, _f) / u;
    }


    // present value of instrument having cash flow c[i] at time u[i]
    template<class T, class F>
    inline F present_value(size_t m, const T* u, const F* c, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        F p{ 0 };

        for (size_t i = 0; i < m; ++i) {
            p += c[i] * pwflat::discount(u[i], n, t, f, _f);
        }

        return p;
    }

    // derivative of present value wrt parallel shift of forward curve
    template<class T, class F>
    inline F duration(size_t m, const T* u, const F* c, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        F d{ 0 };

        for (size_t i = 0; i < m; ++i) {
            d -= u[i] * c[i] * pwflat::discount(u[i], n, t, f, _f);
        }

        return d;
    }

    // derivative of present value wrt parallel shift of forward curve after last curve time
    template<class T, class F>
    inline F partial_duration(size_t m, const T* u, const F* c, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<X>::quiet_NaN()) noexcept
    {
        F d{ 0 };

        // first cash flow past end of forward curve
        size_t i0 = (n == 0) ? 0 : std::lower_bound(u, u + m, t[n - 1]) - u;
        T t0 = (n == 0) ? 0 : t[n - 1];
        for (size_t i = i0; i < m; ++i) {
            d -= (u[i] - t0)*c[i] * pwflat::discount(u[i], n, t, f, _f);
        }

        return d;
    }

} // fms::pwflat
