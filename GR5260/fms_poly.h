// fms_poly.h - Various polynomials
#pragma once

namespace fms::poly {

    // The n-th Bell polynomial B_n(kappa[0],...,kappa[n-1])
    template<class X = double, class K = double>
    inline auto Bell(size_t n, const K* kappa)
    {
        return n + *kappa; //??? return B_n(kappa_0,...,kappa_{n-1});
    }

    // The n-th Hermite polynomial He_n(x)
    template<class X = double>
    inline auto Hermite(size_t n, X x)
    {
        return n + x; //??? Return n-th Hermite polynomial at x.
    }

}