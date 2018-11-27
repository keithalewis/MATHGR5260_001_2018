// correlation.h - correlation matrices
#pragma once
#include <algorithm>
#include <vector>

/*
An n x n correlation matrix, rho, is determined by n unit vectors (e_j).
The i,j entry in the matrix rho_{j,k} = e_j . e_k, the dot product
of the two unit vectors. Conversely, given any correlation matrix
the Cholesky decomposition has row vectors of unit norm.

If the unit vectors lie in a d-dimensional sphere the Cholesky factor has the form

[ 1    0    0 ... 0 ]
[ e_21 e_22 0 ... 0 ]
[ ...
[ e_d1 e_d2 ... e_dd]
[ ...
[ e_n1 e_n2 ... e_nd]

assuming 1-based indexing.

*/

namespace fms {

    template<class X = double>
    class correlation {
        // Cholesky decomposition
        // lower-triangular matrix of unit vectors rows
        std::vector<std::vector<X>> e_;
    public:
        enum layout {
            lower,  // e_00, 0, ...; e_10, e_11, 0 ...;
            packed, // e_10; e_20, e_21; ...
        };
        correlation()
        { }
        correlation(size_t n, size_t d, const X* e, layout type = packed)
            : e_(n)
        {
            if (n == 0)
                return;

            e_[0] = std::vector<X>(d);
            e_[0][0] = X(1);

            size_t off = 0;
            for (size_t i = 1; i < d; ++i) {
                e_[i] = std::vector<X>(d);
                X e2 = X(0);
                for (size_t j = 0; j < i; ++j) {
                    X eij = e[off + j];
                    e_[i][j] = eij;
                    e2 += eij * eij;
                }
                e_[i][i] = sqrt(1 - e2); // NaN if e2 > 1
                off += type == packed ? i : d - 1;
            }
            for (size_t i = d; i < n; ++i) {
                e_[i] = std::vector<X>(d);
                X e2 = X(0);
                for (size_t j = 0; j < d - 1; ++j) {
                    X eij = e[off + j];
                    e_[i][j] = eij;
                    e2 += eij * eij;
                }
                e_[i][d - 1] = sqrt(1 - e2); // NaN if e2 > 1
                off += d - 1;
            }
        }

        // Size of correlation matrix.
        size_t size() const
        {
            return e_.size();
        }

        // Dimension of sphere.
        size_t dimension() const
        {
            return e_[0].size();
        }

        // i,j entry
        X operator()(size_t i, size_t j) const
        {
            return i >= j ? e_[i][j] : X(0);
        }

        // correlation
        X rho(size_t i, size_t j) const
        {
            return std::inner_product(e_[i].begin(), e_[i].end(), e_[j].begin(), X(0));
        }
    };
}

