// fms_njr.h - Normal Jarrow-Rudd model.
#pragma once

namespace fms::prob {

    // Normal Jarrow-Rudd model.
    // Psi(x) = Phi(x) - phi(x) [sum_{n>=3} B_n(0,0,kappa_3,...,kappa_n) H_{n-1}(x)/n!]
    // Where Phi and phi are the standard normal cdf and pdf.
    template<class X = double, class K = double>
    inline X njr_cdf(size_t n, const K* kappa, X x) {
        return n + *kappa + x; //???Implement
    }

}