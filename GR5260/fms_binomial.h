// fms_binomial.h - Binomial measure
/*
The sample space is all sequences of 0's and 1's.
The set space<N> has 2^N elements of distinct intial segments of length N.
The set atom<N> represents that atoms in space <N>.
*/
#pragma once

namespace fms::binomial {

    // {0,1}^N
    template<size_t N>
    struct space {
    };


    template<size_t N>
    struct iterator {
        int k;
        bool operator==(const iterator& i) const
        {
            return k == i.k;
        }
        bool operator!=(const iterator& i) const
        {
            return !operator==(i);
        }
        static iterator begin() { return iterator{0}; }
        static iterator end() { return iterator{N+1}; }
        size_t operator*() const
        {
            return k;
        }
        iterator& operator++() {
            ++k;
            return *this;
        }
    };

    // choose(N,k)/2^N = N/2k * N-1/2(k-1) * ... * N - k + 1/2 / 2^(N-k)
    template<size_t N, class X = double>
    constexpr X probability(size_t k) {
        if (2*k > N)
            return probability(N - k);

        if (k == 0)
            return ldexp(X(1), -N); // 1/2^N

        return (N*probability<N-1,X>(k-1))/(2*k);
    }

}
