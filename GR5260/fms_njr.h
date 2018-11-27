// fms_njr.h - Normal Jarrow-Rudd model.
#pragma once
#include "fms_poly.h"
#include "fms_prob.h"

namespace fms::prob {

    // Normal Jarrow-Rudd model.
    // psi(x) = phi(x) [1 + sum_{n>=3} B_n(0,0,kappa_3,...,kappa_n) H_{n}(x)/n!]
    template<class X = double, class K = double>
    inline X njr_pdf(size_t n, const K* kappa, X x) {
        X psi = normal_pdf(x);
        X sum = X(0);
        X n_ = X(2); // 2!

        for (size_t i = 3; i < n; ++i) {
            n_ *= i;
            sum += poly::Bell(i,kappa)*poly::Hermite(i,x)/n_;
        }

        return psi*(1 + sum);
    }
    // Psi(x) = Phi(x) - phi(x) [sum_{n>=3} B_n(0,0,kappa_3,...,kappa_n) H_{n-1}(x)/n!]
    // Where Phi and phi are the standard normal cdf and pdf.
    template<class X = double, class K = double>
    inline X njr_cdf(size_t n, const K* kappa, X x) {
        X Psi = normal_cdf(x);
        X sum = X(0);
        X n_ = X(2); // 2!

        for (size_t i = 3; i < n; ++i) {
            n_ *= i;
            sum += poly::Bell(i, kappa)*poly::Hermite(i-1,x)/n_;
        }

        return Psi - normal_pdf(x)*sum;
    }

}