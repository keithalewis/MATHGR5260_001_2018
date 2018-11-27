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
            X sqrdt = sqrt(u - t);
            std::normal_distribution<X> Z;

            // B += e . dB
            // [ e_00 0 ... 0 ] [dB_1]
            // [ e_10 e_11 ...] [...]
            // [ ....           [dB_d
            // B[j] += e_j0 dB_0 + ... e_jd dB_d
            for (size_t k = 0; k < e.dimension(); ++k) {
                auto dB = sqrdt*Z(r);
                for (size_t j = 0; j < B.size(); ++j) {
                    B[j] += e(j, k)*dB;
                }
            }

            t = u;
        }
        // Size of Brownian sample vector.
        size_t size() const
        {
            return e.size();
        }
        // Dimension of correlation.
        size_t dimension() const
        {
            return e.dimension();
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
