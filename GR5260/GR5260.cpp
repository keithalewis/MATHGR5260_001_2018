// GR5260.cpp - test program
#include <cassert>
#include <algorithm>
#include <chrono>
#include <functional>
#include "fms_analytic.h"
#include "fms_black.h"
#include "fms_njr.h"
#include "fms_poly.h"
#include "fms_root1d_newton.h"
#include "fms_bootstrap.h"

using namespace fms;

inline auto timer(const std::function<void(void)>& f, size_t count = 1)
{
    auto start = std::chrono::system_clock::now();
    while (count--) {
        f();
    }
    auto stop = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = stop-start;

    return elapsed.count(); // duration in seconds
}

template<class X>
void test_fms_poly_Hermite()
{
    using fms::poly::Hermite;
    std::vector<X> x;

    for (auto x_ = X(-2); x_ <= X(2); x_ += X(0.1)) {
        x.push_back(x_);
    }
    
    for (auto x_ : x) {
        assert(1 == Hermite(0, x_));
        assert(x_ == Hermite(1, x_));
        assert(x_*x_ - 1 == Hermite(2, x_));
        assert(x_*(x_*x_ - 1) - 2 * x_ == Hermite(3, x_));
        //??? put in tests for n = 4 and n = 5
    }
}

template<class X>
void test_fms_poly_Bell()
{
    using fms::poly::Bell;
    X kappa[] = { X(1), X(2), X(3), X(4) };

    assert(1 == Bell(0, kappa));
    X B1 = kappa[0];
    assert(B1 == Bell(1, kappa));
    X B2 = B1 * kappa[0] + kappa[1];
    assert(B2 == Bell(2, kappa));
    X eps;
    X B3 = B2*kappa[0] + 2*B1*kappa[1] + kappa[2]; 
    eps = B3 - Bell(3,kappa);
    assert (B3 == Bell(3,kappa));
    X B4 = B3*kappa[0] + 3*B2*kappa[1] + 3*B1*kappa[2] + kappa[3];
    eps = B4 - Bell(4,kappa);
    assert (B4 == Bell(4,kappa));

    for (size_t i = 0; i < 5; ++i) {
        auto Bi = Bell(i,kappa);
        assert (Bi == fms::poly::Bell2(i,kappa));
        assert (Bi == fms::poly::Bell3(i,kappa));
    }
    double secs;
    secs = timer([&kappa]() { Bell(4,kappa); }, 10000);
    secs = secs;
    secs = timer([&kappa]() { fms::poly::Bell2(4,kappa); }, 10000);
    secs = secs;
    secs = timer([&kappa]() { fms::poly::Bell3(4,kappa); }, 10000);
    secs = secs;
}

template<class X>
void test_fms_prob_njr()
{
    X x = 0;
    std::vector<X> kappa(4);

    for (int i = 0; i < 1000; ++i) {
        fms::prob::njr_cdf(4,kappa.data(),x + i/10000.);
    }
    x = x;

    // Independent check.
    // X = Normal(0, 1) + Poisson(lambda)
    // P(X <= x) = sum_{k>=0} P(X + k <= x) exp(-lambda)lambda^k/k!
}

template<class X>
void test_fms_root1d_newton()
{
    size_t n = 0;
    X a = X(5);
    std::function<X(X)> f = [a](X x) { return x * x - a; };
    std::function<X(X)> df = [](X x) { return 2 * x; };
    X x0 = X(2);
    root1d::newton_solver<X,X> ns(X(x0), f, df);

    // by hand to see how many iterations
    X x_;
    do {
        ++n;
        x_ = ns.next();
    } while (!ns.done());

    // should get same solution
    X x = ns.solve();
    assert(x == x_);

    bool thrown = false;
    try {
        root1d::newton_solver<X, X, 2>(X(x0), f, df).solve();
    }
    catch (...) {
        thrown = true;
    }
    assert(thrown);
}

template<class X>
std::vector<X> sequence(X start, X stop, X step = X(1))
{
    std::vector<X> v;

    while (start < stop) {
        v.push_back(start);
        start += step;
    }

    return v;
}

template<class X>
void test_fms_black_value() 
{
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

    int n = 15; // 10
    assert (fabs(dv - dv0) < n*eps);
}
template<class X>
void test_fms_black_vega()
{
    X f = X(100), sigma = X(.2), k = X(100), t = X(0.25);
    X eps = sqrt(std::numeric_limits<X>::epsilon());

    auto _v = black::value(f, sigma - eps, k, t);
    auto v_ = black::value(f, sigma + eps, k, t);

    auto dv = (v_ - _v) / (2 * eps);
    auto dv0 = black::vega(f, sigma, k, t);

    auto diff = dv - dv0;
    auto nn = diff / eps;

    int n = 5; // 3
    assert(fabs(dv - dv0) < n*eps);
}

// multiple of machine epsilon
template<class X>
struct implied {
    static const int n;
};
template<>
struct implied<double> {
    static const int n = 3; 
};
template<>
struct implied<float> {
    static const int n = 4; 
};

template<class X>
void test_fms_black_implied()
{
    std::vector<X> sigma;
    // vol from 1% to 100% in 10% increments
    for (X x = X(0.1); x <= X(1); x += X(0.1)) {
        sigma.push_back(x);
    }
    X f = X(100), k = X(100), t = X(0.25);
    X eps = std::numeric_limits<X>::epsilon();
    
    for (auto sigmai : sigma) {
        X v = black::value(f, sigmai, k, t);
        X s = black::implied(f, v, k, t);
        X eps_ = s - sigmai;
        assert(fabs(eps_) <= implied<X>::n*eps);
    }
}


template<class X>
void test_fms_black() 
{
    test_fms_black_value<X>();
    test_fms_black_delta<X>();
    test_fms_black_vega<X>();
    test_fms_black_implied<X>();
}

template<class X>
void test_fms_analytic()
{
    using fms::analytic;

    analytic<X> x1(1);
    assert (x1.size() == 1);
}

int main()
{
    test_fms_analytic<double>();

    test_fms_poly_Hermite<double>();
    test_fms_poly_Hermite<float>();

    test_fms_poly_Bell<double>();
    test_fms_poly_Bell<float>();

    test_fms_prob_njr<double>();

    test_fms_root1d_newton<double>();
    test_fms_root1d_newton<float>();

    test_fms_black<double>();
    test_fms_black<float>();
}