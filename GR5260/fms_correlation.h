// correlation.h - correlation matrices
#pragma once
#include <vector>

/*
An n x n correlation matrix, rho, is determined by n unit vectors (e_j).
The i,j entry in the matrix rho_{j,k} = e_j . e_k, the dot product
of the two unit vectors. Conversely, given any correlation matrix
the Cholesky decomposition has row vectors of unit norm.
*/

namespace fms {

    template<class X = double>
    class correlation {
        // Cholesky decomposition
        // lower-triangular matrix of unit vectors rows
        std::vector<std::vector<X>> e_;
    public:
        enum layout {
            lower,  // e_10, 0, ...; e_20, e_21, 0 ...;
            packed, // e_10; e_20, e_21; ...
        };
        correlation()
        { }
        correlation(size_t n, const X* e, layout type = packed)
            : e_(n)
        {
            if (n == 0)
                return;

            e_[0] = std::vector<X>{X(1)};

            for (size_t i = 1; i < n; ++i) {
                e_[i] = std::vector<X>(i+1);
                X e2 = X(0);
                for (size_t j = 0; j < i; ++j) {
                    size_t off = (type != packed ? (i - 1)*(n - 1) : (i*(i-1))/2);
                    X eij = e[off + j];
                    e_[i][j] = eij;
                    e2 += eij * eij;
                }
                e_[i][i] = sqrt(1 - e2); // NaN if e2 > 1
            }
        }

        size_t size() const
        {
            return e_.size();
        }

        // i,j entry
        X operator()(size_t i, size_t j) const
        {
            return i >= j ? e_[i][j] : X(0);
        }

        // correlation
        X rho(size_t i, size_t j) const
        {
            //ensure(i < n && j < n);

            if (i > j) {
                std::swap(i, j);
            }

            X r{0};
            const auto& e_i = e_[i];
            const auto& e_j = e_[j];
            for (size_t k = 0; k <= i; ++k) {
                r += e_i[k]*e_j[k];
            }

            return r;
        }
    };
}

