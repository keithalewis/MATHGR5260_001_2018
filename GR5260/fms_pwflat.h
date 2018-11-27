// fms_pwflat.h - piecewise flat curve
/*
    f(t) = f[i] if t[i-1] < t <= t[i]
         = _f   if t > t[n-1]
    and undefined if t < 0
    Note f(0) = f[0] and f(t[i]) = f[i] for all i

    |                                   _f
    |        f[1]             f[n-1] (--------
    | f[0] (----- ...       (--------]
    [------]      ... ------]
    |
    0----t[0]--- ...  ---t[n-2]---t[n-1]
*/
#pragma once
#include <cmath>     // exp
#include <algorithm> // adjacent_find
#include <limits>    // quiet_Nan()
#include <numeric>   // upper/lower_bound
#include "fms_fixed_income_instrument.h"

namespace fms::pwflat {

    // strictly increasing values
    template<class I>
    inline bool strictly_increasing(I b, I e)
    {
        return e == std::adjacent_find(b, e, std::greater_equal{});
    }
    template<class T>
    inline bool strictly_increasing(size_t n, const T* t) noexcept
    {
        return strictly_increasing(t, t + n);
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
    template<class U, class C, class T, class F>
    inline F present_value(size_t m, const U* u, const C* c, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        F p{ 0 };

        for (size_t i = 0; i < m; ++i) {
            p += c[i] * pwflat::discount<T,F>(u[i], n, t, f, _f);
        }

        return p;
    }

    // derivative of present value wrt parallel shift of forward curve
    template<class U, class C, class T, class F>
    inline F duration(size_t m, const U* u, const C* c, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        F d{ 0 };

        for (size_t i = 0; i < m; ++i) {
            d -= u[i] * c[i] * pwflat::discount(u[i], n, t, f, _f);
        }

        return d;
    }

    // derivative of present value wrt parallel shift of forward curve after last curve time
    template<class U, class C, class T, class F>
    inline F partial_duration(size_t m, const U* u, const C* c, size_t n, const T* t, const F* f, 
        const F& _f = std::numeric_limits<F>::quiet_NaN()) noexcept
    {
        F d{ 0 };

        // first cash flow past end of forward curve
        size_t i0 = (n == 0) ? 0 : std::lower_bound(u, u + m, t[n - 1]) - u;
        T t0 = (n == 0) ? 0 : t[n - 1];
        for (size_t i = i0; i < m; ++i) {
            d -= (u[i] - t0)*c[i] * pwflat::discount<T,F>(u[i], n, t, f, _f);
        }

        return d;
    }

    // NVI curve interface.
    template<class T = double, class F = double>
    class curve {
        size_t n;
        T* t;
        F* f;
        F _f;
    public:
        typedef T time_type;
        typedef F rate_type;
        // constant value
        curve(F f)
            : curve(0, nullptr, nullptr, f)
        { }
        curve(size_t n = 0, T* t = 0, F* f = 0, F _f = std::numeric_limits<F>::quiet_NaN())
            : n(n), t(t), f(f), _f(_f)
        { }
        virtual ~curve() {}

        // convenience functions
        bool operator==(const curve& i) const
        {
            return size() == i.size()
                && std::equal(time(), time() + size(), i.time())
                && std::equal(rate(), rate() + size(), i.rate());
        }
        bool operator!=(const curve& i) const
        {
            return !operator==(i);
        }

        size_t   size() const { return _size(); }
        const T* time() const { return _time(); }
        const F* rate() const { return _rate(); }

        F value(T u) const
        {
            return pwflat::value<T,F>(u, size(), time(), rate(), _f);
        }
        F operator()(T u) const
        {
            return value(u);
        }
        F discount(T u) const
        {
            return pwflat::discount<T,F>(u, size(), time(), rate(), _f);
        }
        F spot(T u) const
        {
            return spot<T,F>(u, size(), time(), rate(), _f);
        }
        F forward(T u) const
        {
            return pwflat::forward<T,F>(u, size(), time(), rate(), _f);
        }
        template<class U, class C>
        F present_value(const fixed_income::instrument<U,C>& i)
        {
            return pwflat::present_value<U,C,T,F>(i.size(), i.time(), i.cash(), size(), time(), rate(), _f);
        }
        template<class U, class C>
        F duration(const fixed_income::instrument<U,C>& i)
        {
            return pwflat::duration<U,C,T,F>(i.size(), i.time(), i.cash(), size(), time(), rate(), _f);
        }
    private:
        virtual size_t _size() const
        {
            return n;
        }
        virtual T* _time() const
        {
            return t;
        }
        virtual F* _rate() const
        {
            return f;
        }
    };

} // fms::pwflat
