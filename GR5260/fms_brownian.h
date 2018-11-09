// fms_brownian.h - Brownian motion
#pragma once
#include <algorithm>
#include <random>
#include <vector>
#include "fms_correlation.h"

namespace fms {

    // n-dimensional correlated Brownian samples
    template<class T = double, class X = double>
    class brownian {
        T t;
        fms::correlation<X> e;
        std::vector<X> B;
    public:
        brownian(const fms::correlation<X>& e)
            : t(T(0)), e(e), B(e.size())
        {
            reset();
        }

        void reset()
        {
            t = T(0);
            std::fill(B.begin(), B.end(), (X(0)));
        }
        // Generate B_u using a random engine
        template<class R> // random engine
        void advance(X u, R& r)
        {
            std::normal_distribution<X> Z(0, sqrt(u - t));

            // B += e . dB
            for (size_t i = 0; i < B.size(); ++i) {
                auto dB = Z(r);
                for (size_t j = i; j < B.size(); ++j) {
                    B[j] += e(j, i)*dB;
                }
            }

            t = u;
        }
        // number of dimensions
        size_t size() const
        {
            return B.size();
        }
        const X* data() const
        {
            return B.data();
        }
        X operator[](size_t i) const
        {
            return B[i];
        }
        T time() const
        {
            return t;
        }
    };

}
