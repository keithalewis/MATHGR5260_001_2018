// fms_poly.h - Various polynomials
#pragma once

namespace fms::poly {

    // The Bell polynomial B_n(kappa_1,...,kappa_n) satisfy B_0 = 1 and
    //   B_{n+1}(kappa_1, ..., kappa_{n+1}) = sum_{k=0}^n C(n,k) B_{n-k}(kappa_1,...,kappa_{n-k}) kappa_{k+1}.
    // Replacing n by n-1 and renaming kappa_j to kappa_{j-1} yields
    //   B_n(kappa_0, ..., kappa_{n-1}) = sum_{k=0}^{n-1} C(n-1,k) B_{n-1-k}(kappa_0,...,kappa_{n-1-k}) kappa_k.
    template<class K = double>
    inline auto Bell(size_t n, const K* kappa)
    {
        if (n == 0) {
            return K(1);
        }

        K n_ = K(n - 1);
        K C = 1; // C(n-1,0);
        K B = 0;

        for (size_t k = 0; k < n; ++k) {
            B += C * Bell(n - 1 - k, kappa) * kappa[k];
            C *= n_;
            C /= k + 1; // C(n - 1, k)
            n_ = n_ - 1;
        }

        return B;
    }
    // B_n(kappa_0, ..., kappa_{n-1}) 
    //   = sum_{k=0}^{n-1} C(n-1,k) B_{n-1-k} kappa_k.
    //   = sum_{k=0}^{n-1} C(n-1,n-1-k) B_{n-1-k} kappa_k.
    //   = sum_{k=0}^{n-1} C(n-1,k)B_k kappa_{n-1-k} 
    template<class K = double>
    inline auto Bell2(size_t n, const K* kappa)
    {
        if (n == 0) {
            return K(1);
        }

        K n_ = K(n - 1);
        K C = 1; // C(n-1,0);
        K B = 0;

        for (size_t k = 0; k < n; ++k) {
            B += C * Bell2(k, kappa) * kappa[n-1-k];
            C *= n_;
            C /= k + 1; // C(n - 1, k)
            n_ = n_ - 1;
        }

        return B;
    }
    // Bell with B_1, ..., B_m already computed. B_ should have room for n+1 values.
    template<class K = double>
    inline auto Bell3(size_t n, const K* kappa)
    {
        static const K* okappa = nullptr;
        static std::vector<K> B_;

        if (kappa != okappa) {
            okappa = kappa;
            B_.resize(0);
        }

        if (n < B_.size()) {
            return B_[n]; // already computed
        }

        if (n == 0) {
            B_.push_back(K(1));

            return K(1);
        }
 
        K n_ = K(n - 1);
        K C = 1; // C(n-1,0);
        K B = 0;

        for (size_t k = 0; k < B_.size(); ++k) {
            B += C * B_[k] * kappa[n-1-k];
            C *= n_;
            C /= k + 1; // C(n - 1, k)
            n_ = n_ - 1;
        }
        B_.push_back(B);
        for (size_t k = B_.size(); k < n; ++k) {
            B_.push_back(Bell(k, kappa));
        }

        return B_[n];
    }


    // The Hermite polynomials satisfy He_0(x) = 1, He_1(x) = x, and 
    //   He_{n+1}(x) = x He_n(x) - n He_{n-1}(x). Replacing n by n - 1
    //   He_n(x) = x He_{n-1}(x) - (n-1} He_{n-2}(x).
    template<class X = double>
    inline X Hermite(size_t n, X x)
    {
        if (n == 0) {
            return 1;
        }

        if (n == 1) {
            return x;
        }

        return x * Hermite(n - 1, x) - (n - 1)*Hermite(n - 2, x);
    }

}