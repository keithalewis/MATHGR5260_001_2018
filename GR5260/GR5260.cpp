// GR5260.cpp - test program
#include <cassert>
#include "root1d.h"
#include "fms_black.h"

using namespace fms;

template<class X>
void test_fms_black_value() {
    X f = X(100), sigma = X(.2), k = X(100), t = X(0.25);

    auto v = black::value(f, sigma, k, t);
    assert(fabs(v - 4.0) < .02);

    auto s = sigma * sqrt(t);
    auto v2 = black::value(f, s, k);
    assert(v == v2);
}
template<class X>
void test_fms_black_delta()
{
    X f = X(100), sigma = X(.2), k = X(100), t = X(0.25);
    X eps = sqrt(std::numeric_limits<X>::epsilon());

    auto _v = black::value(f - eps, sigma, k, t);
    auto v_ = black::value(f + eps, sigma, k, t);

    auto dv = (v_ - _v)/(2*eps);
    auto dv0 = black::delta(f, sigma, k, t);

    int n = 0; //??? find the smallest n that make the following test pass
    assert (fabs(dv - dv0) < n*eps);
}

template<class X>
void test_fms_black() {
    test_fms_black_value<X>();
    test_fms_black_delta<X>();
}

int main()
{
    test_fms_black<double>();
    test_fms_black<float>();
}