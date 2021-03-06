// GR5260.cpp - test program
#include <cassert>
#include <algorithm>
#include <chrono>
#include <functional>
#include <random>
#include "fms_analytic.h"
#include "fms_black.h"
#include "fms_brownian.h"
#include "fms_njr.h"
#include "fms_poly.h"
#include "fms_root1d_newton.h"
#include "fms_bootstrap.h"
#include "fms_fixed_income.h"
#include "fms_ho_lee.h"
#include "fms_swaption.h"

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

// [f() + ... f()]/n
template<class X>
inline X mean(const std::function<X()>& f, size_t n)
{
     X m = 0;

    for (size_t i = 1; i <= n; ++i)
        m += (f() - m)/i;

    return m;
}

template<class X>
void test_mean()
{
    for (size_t n = 1; n < 10; ++n) {
        X k = 0;
        std::function<X()> f = [n,&k]() { return ++k/n; };
        X m = mean(f, n);
        // (1/n + ... + n/n)/n = n(n+1)/2n^2
        assert (m == X(n + 1)/(2*n));
    }

    std::default_random_engine dre;
    std::uniform_real_distribution<X> u;

    size_t N = 10000;
    std::function<X()> f = [u,&dre]() { return u(dre); };
    X m = mean(f, N);
    assert (fabs(m - X(0.5)) < X(1)/sqrt(N));
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
        X H2 = x_*x_ - 1;
        assert(H2 == Hermite(2, x_));
        X H3 = x_*H2 - 2*x_;
        assert(H3 == Hermite(3, x_));
        X H4 = x_*H3 - 3*H2;
        assert(H4 == Hermite(4, x_));
        X H5 = x_*H4 - 4*H3;
        assert(H5 == Hermite(5, x_));
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

//    auto diff = dv - dv0;
//    auto nn = diff / eps;

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
    
    {
        analytic<X> x(1);
        assert (x.order() == 1);
        assert (x[0] == X(0));
        assert (x(0) == X(0));
        analytic<X> x2(x);
        assert (x == x2);
        x = x2;
        assert (x == x2);
        x = X(3);
        assert (x[0] == X(3));
        x += x;
        assert (x[0] == X(6));
    }
    {
        analytic<X> x(2);
        assert (x.order() == 2);
        assert (x[0] == X(0));
        assert (x[1] == X(0));
        assert (x(0) == X(0));
        assert (x(1) == X(0));
        analytic<X> x2(x);
        assert (x == x2);
        x = x2;
        assert (x == x2);
    }
    {
        analytic<X> x{X(2),X(1)};
        assert (x.order() == 2);
        assert (x[0] == X(2));
        assert (x[1] == X(1));
        assert (x(0) == X(2));
        assert (x(1) == X(1));
        x += x;
        assert (x[0] == X(4));
        assert (x[1] == X(2));
        assert (x(0) == X(4));
        assert (x(1) == X(2));

        analytic x2{X(2)};
        assert (x2 != x);
        x += x2;
        assert (x[0] == X(6));
        assert (x[1] == X(2));
    }
    {
        analytic x{X(2),X(1),X(0)};
        assert (x[0] == X(2));
        assert (x[1] == X(1));
        assert (x[2] == X(0));
        
        auto x2 = x;
        x2 *= x;
        // (2 + J)*(2 + J) = 4 + 2*2J + J^2
        assert (x2[0] == X(4));
        assert (x2[1] == X(4));
        assert (x2[2] == X(1));

        auto x3 = x*x;
        assert (x3 == x2);

        analytic x4{X(2),X(1)};
        auto x5 = x*x4;
        assert (x5 == x2);
    }
}

template<class X>
void test_fms_pwflat()
{
    using namespace fms::pwflat;

    std::vector<X> t{1.,2.,3.}, f{.1,.2,.3};
    std::vector<X> t_2{ 1. }, f_2{ .1 };

    { // monotonic
        ensure (strictly_increasing(std::begin(t), std::end(t)));
        ensure (strictly_increasing(std::begin(f), std::end(f)));
        X f2 = f[2];
        f[2] = -1;
        ensure (!strictly_increasing(std::begin(f), std::end(f)));
        f[2] = f2;
        ensure (!strictly_increasing(std::rbegin(f), std::rend(f)));
    }
    { // forward
      //0, 0, null, null, null
        ensure (isnan(value<int,X>(0, 0, nullptr, nullptr)));
        //1, 0, null, null, null
        ensure(isnan(value<int, X>(1, 0, nullptr, nullptr)));
        //-1, 0, null, null, null
        ensure(isnan(value<int, X>(-1, 0, nullptr, nullptr)));
        //-1, 0, null, null, 0.2
        ensure(isnan(value<int, X>(-1, 0, nullptr, nullptr, 0.2)));

        int u;
        u = 1;
        X x{ 0.2 }, x_;
        //1, 0, null, null, 0.2
        x_ = fms::pwflat::value<int, X>(u, 0, nullptr, nullptr, x);
        ensure(x_ == x);

        X u_ [] = { -1, 0, 0.5, 1, 1.5 };
        X a_ [] = { 0, 0.1, 0.1, 0.1, 0.2 };

        for (int i = 0; i < 5; i++) {
            if (i == 0 || i == 4) {
                ensure(isnan(value<X, X>(u_[i], t_2.size(), t_2.data(), f_2.data())));
            }
            else {
                x_ = fms::pwflat::value<X, X>(u_[i], t_2.size(), t_2.data(), f_2.data());
                ensure(x_ == a_[i]);
            }
        }

        for (int i = 0; i < 5; i++) {
            if (i == 0) {
                ensure(isnan(value<X, X>(u_[i], t_2.size(), t_2.data(), f_2.data(), 0.2)));
            }
            else {
                x_ = fms::pwflat::value<X, X>(u_[i], t_2.size(), t_2.data(), f_2.data(), 0.2);
                ensure(x_ == a_[i]);
            }
        }

        for (int i = 0; i < 3; ++i)
            ensure (f[i] == value(t[i], t.size(), t.data(), f.data()));
    }
    { // integral
        X u;
        u = -1;
        ensure (isnan(integral(u, t.size(), t.data(), f.data())));
        u = 4;
        ensure (isnan(integral(u, t.size(), t.data(), f.data())));
        u = 0;
        ensure (0 == integral(u, t.size(), t.data(), f.data()));
        u = 0.5;
        ensure (.1*.5 == integral(u, t.size(), t.data(), f.data()));
        u = 1;
        ensure (.1 == integral(u, t.size(), t.data(), f.data()));
        u = 1.5;
        ensure (.1 + .2*.5 == integral(u, t.size(), t.data(), f.data()));
        u = 2.5;
        ensure (.1 + .2 + .3*.5 == integral(u, t.size(), t.data(), f.data()));
        u = 3;
        ensure (fabs(.1 + .2 + .3 - integral(u, t.size(), t.data(), f.data())) < 1e-10);
        //		ensure (.1 + .2 + .3 != .6); 
    }
    { // discount
        X u_[] = { -.5, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5 };
        X f_[] = {0, 0, .05, .1, .2, .3, .45, .6, .7};
        for (int i = 0; i < 9; i++) {
            if (i == 0 || i == 8) {
                ensure(isnan(discount(u_[i], t.size(), t.data(), f.data())));
            } 
            else {
                ensure(fabs(exp(-f_[i]) - discount(u_[i], t.size(), t.data(), f.data())) < 1e-10);
            }
        }

        for (int i = 0; i < 9; i++) {
            if (i == 0) {
                ensure(_isnan(discount(u_[i], t.size(), t.data(), f.data(), 0.2)));
            }
            else {
                ensure(fabs(exp(-f_[i]) - discount(u_[i], t.size(), t.data(), f.data(), 0.2)) < 1e-10);
            }
        }
    }
    { // spot
        X u_[] = { -.5, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5 };
        X f_[] = { .1, .1, .1, .1, .2/1.5, .3/2, .45/2.5, .6/3, .7/3.5 };
        for (int i = 0; i < 9; i++) {
            if (i == 8) {
                ensure(isnan(spot(u_[i], t.size(), t.data(), f.data())));
            }
            else {
                ensure(fabs(f_[i] - spot(u_[i], t.size(), t.data(), f.data())) < 1e-10);
            }
        }

        for (int i = 0; i < 9; i++) {
            ensure(fabs(f_[i] - spot(u_[i], t.size(), t.data(), f.data(), 0.2)) < 1e-10);
        }
    }
    { // present_value
        X u_[] = { 0, 1, 2, 3, 4};
        X d_[] = { 0,
            discount(u_[1], t.size(), t.data(), f.data(), 0.2),
            discount(u_[2], t.size(), t.data(), f.data(), 0.2),
            discount(u_[3], t.size(), t.data(), f.data(), 0.2),
            discount(u_[4], t.size(), t.data(), f.data(), 0.2)
        };
        X c_[] = { 0, 1, 2, 3, 4 };

        //ensure(isnan(present_value(1, u_, c_, t.size(), t.data(), f.data())));
        //ensure(isnan(present_value(1, u_, c_, t.size(), t.data(), f.data(), 0.2)));

        X sum = 0;
        for (int i = 0; i < 5; i++) {
            sum += c_[i] * d_[i];
            if (i == 4) {
                X tmp = present_value<X, X>(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2);
                ensure(tmp == tmp);
                ensure(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2)) < 1e-10);
                ensure(isnan(present_value(i + 1, u_, c_, t.size(), t.data(), f.data())));
            }
            else {
                X tmp = present_value<X, X>(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2);
                ensure(tmp == tmp);
                ensure(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data(), 0.2)) < 1e-10);
                ensure(fabs(sum - present_value(i + 1, u_, c_, t.size(), t.data(), f.data())) < 1e-10);
            }
        }

    }
}

template<class X>
void test_fms_fixed_income_zero()
{
    using fms::fixed_income::zero;
    using fms::pwflat::bootstrap;
    using fms::pwflat::curve;

    zero z(X(1));
    zero z2(z); // copy constructor
    assert (z == z2);
    z = z2; // copy assignment
    assert (z == z2);
    assert (z.size() == 1);
    auto u = z.time();
    assert (u[0] == 1);
    auto c = z.cash();
    assert (c[0] == 1);

    X r = X(0.01);
    std::vector<X> T, F;

    for (X u0 = 1; u0 < 10; u0 = u0 + X(1)) {
        auto [t0,f0] = bootstrap(exp(-r*u0), zero(u0), curve(T.size(), T.data(), F.data()));
        T.push_back(t0);
        F.push_back(f0);
    }

    curve f(T.size(), T.data(), F.data());
    X flo = std::numeric_limits<X>::max();
    X fhi = -std::numeric_limits<X>::max();
    for (X u0 = 0; u0 < 10; u0 = u0 + X(0.1)) {
        X eps;
        eps = f(u0) - r;
        if (eps > fhi) fhi=eps;
        if (eps < flo) flo=eps;
    }
    assert (flo > -std::numeric_limits<X>::epsilon());
    assert (fhi < std::numeric_limits<X>::epsilon());
}

template<class X>
void test_fms_pwflat_bootstrap()
{
    using fms::fixed_income::frequency;
    using fms::fixed_income::zero;
    using fms::fixed_income::cash_deposit;
    using fms::fixed_income::forward_rate_agreement;
    using fms::fixed_income::interest_rate_swap;
    using fms::pwflat::bootstrap;
    using fms::pwflat::curve;

    X r = 0.02;
    auto f0 = curve<X,X>(0, nullptr, nullptr, r);
    assert (f0(100) == r);

    std::vector<X> T, F;

    auto cd0 = cash_deposit<X,X>(0.25, r);
    auto cd1 = cash_deposit<X,X>(0.5, r);
    auto cd2 = cash_deposit<X,X>(1.0, r);

    auto fra0 = forward_rate_agreement<X,X>(1.,1.25,r);
    auto fra1 = forward_rate_agreement<X,X>(1.25,1.5,r);
    auto fra2 = forward_rate_agreement<X,X>(1.5,1.75,r);
    auto fra3 = forward_rate_agreement<X,X>(1.75,2.,r);

    auto irs0 = interest_rate_swap<X,X>(3, r, frequency::semiannual);
    auto irs1 = interest_rate_swap<X,X>(5, r, frequency::quarterly);
    auto irs2 = interest_rate_swap<X,X>(10, r, frequency::monthly);

    {    
       auto [t,f] = bootstrap(f0.present_value(cd0), cd0, curve(T.size(), T.data(), F.data()));
       T.push_back(t);
       F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(cd1), cd1, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(cd2), cd2, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(fra0), fra0, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(fra1), fra1, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(fra2), fra2, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(fra3), fra3, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(irs0), irs0, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(irs1), irs1, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }
    {    
        auto [t,f] = bootstrap(f0.present_value(irs2), irs2, curve(T.size(), T.data(), F.data()));
        T.push_back(t);
        F.push_back(f);
    }

    auto f = curve<X,X>(T.size(), T.data(), F.data());
    X flo = std::numeric_limits<X>::max();
    X fhi = -std::numeric_limits<X>::max();
    for (X u0 = 0; u0 < 10; u0 = u0 + X(0.1)) {
        X eps;
        eps = f(u0) - r;
        if (eps > fhi) fhi=eps;
        if (eps < flo) flo=eps;
    }
    assert (flo > -3*std::numeric_limits<X>::epsilon());
    assert (fhi < 3*std::numeric_limits<X>::epsilon());
}

template<class X>
void test_fms_brownian()
{
    // Show corr(B_1[j], B_1[k]) = rho_{j,k}
    X e[] = {X(0.1), X(0.2), X(0.3)};
    fms::correlation<X> corr(3, 3, e); // packed
    // Cholesky decompostion is
    // [ 1   0               0                      ]
    // [ 0.1 sqrt(1 - 0.1^2) 0                      ]
    // [ 0.2 0.3             sqrt(1 - 0.2^2 - 0.3^2)]
    fms::brownian<X> B(corr);
    std::default_random_engine dre;

    size_t N = 10'000; // number of simulations
    for (size_t j = 0; j < B.size(); ++j) {
        for (size_t k = 0; k < B.size(); ++k) {
            X rho = corr.rho(j, k);
            // corr(B_1[j],B_1[k]) = Cov(B_1[j],B_1[k]) = E B_1[j] B_1[k]
            std::function<X()> f = [j,k,&B,&dre]() { 
                B.reset(); 
                B.advance(1,dre); 
                return B[j]*B[k]; 
            };
            assert (fabs(mean(f,N) - rho) < X(3)/sqrt(N));
        }
    }
}
/*
// int_0^1 B_s ds
template<class X, class R>
inline X intB(size_t n, R& rng)
{
    X dt = X(1)/n;
    std::normal_distribution<X> dB(0,sqrt(dt));

    X B = X(0);
    X intBds = X(0);
    for (size_t i = 0; i < n; ++i) {
        B += dB(rng);
        intBds += B*dt;
    }

    return intBds;
}

// Var int_0^1 B_s ds = 1/3
template<class X>
inline void test_intB()
{
    std::default_random_engine dre;
    unsigned int now = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
    dre.seed(now);
    size_t n = 100;
    size_t N = 5000;
    std::function<X()> f = [n,N,&dre]() { X B = intB<X>(n,dre); return B*B; };
    auto var = mean(f, N);
    assert (fabs(var - X(1)/3) < 2/sqrt(N));
}
*/
template<class X>
void test_ElogD_()
{
    // log D_t(u) =  -sigma(u - t)B_t - int_t^u [phi(s) - sigma^2(u - s)^2/2] ds).
}

template<class X>
void test_fms_ho_lee()
{
    X u = X(1), v = X(1.25);
    X f = X(0.02), k = X(0.02), sigma = X(0.20);
    X dcf = v - u;
    pwflat::curve<X,X> F(f); // constant curve
    X Du = F.discount(u);
    X Dv = F.discount(v);

    X p = ho_lee::floor(k, dcf, u, v, Du, Dv, sigma);

    X eld = ho_lee::ElogD_(u, v, Du, Dv, sigma);
    X vld = ho_lee::VarlogD_(u, v, sigma);
    std::default_random_engine dre;
    std::normal_distribution<double> Z(eld, sqrt(vld));

    std::function<X()> p_ = [k,dcf,&Z,&dre]() {
        auto r = (1/exp(Z(dre)) - 1)/dcf; // Ho-Lee forward rate
        return std::max<X>(k - r, 0);
    };
    size_t N = 10000;
    X ep = mean<X>(p_, N);
    // should be within 2 standard deviations
    //assert(fabs(p - ep)/p < 2/sqrt(X(N)));
}

template<class X>
#pragma warning(disable: 4456) // declaration of 'corr' hides previous local declaration)
void test_fms_correlation()
{
    double eps = std::numeric_limits<double>::epsilon();
    //using fms::correlation<double>::packed();
    {
        fms::correlation<> corr;
    }
    {
        fms::correlation<> corr(1, 1, nullptr);
        ensure (corr.size() == 1);
        ensure (corr.rho(0, 0) == 1);
    }
    {
        double rho = 0.5;
        fms::correlation<> corr(2, 2, &rho, fms::correlation<>::packed);
        ensure (corr.size() == 2);
        ensure (corr.rho(0, 0) == 1);
        ensure (corr.rho(0, 1) == 0.5);
        ensure (corr.rho(1, 0) == 0.5);
        ensure (fabs(corr.rho(1, 1) - 1) <= eps);
    }
    {
        double rho[] = {0.5,  0.4, 0.3};
        fms::correlation<> corr(3, 3, rho);
        ensure (corr.size() == 3);
        ensure (corr.rho(0, 0) == 1);
        ensure (corr.rho(0, 1) == 0.5);
        ensure (corr.rho(0, 2) == 0.4);
        ensure (corr.rho(1, 0) == 0.5);
        ensure (fabs(corr.rho(1, 1) - 1) <= eps);
        ensure (corr.rho(1, 2) == 0.5*0.4 + sqrt(1-0.5*0.5)*0.3);
        ensure (fabs(corr.rho(2, 0) - 0.4) <= eps);
        ensure (fabs(corr.rho(2, 1) - (0.5*0.4 + sqrt(1-0.5*0.5)*0.3)) <= eps);
        ensure (fabs(corr.rho(2, 2) - 1) <= eps);
    }
    {
        double rho[] = {0.5, 0, 0.4, 0.3};
        fms::correlation<> corr(3, 3, rho, fms::correlation<>::lower);
        ensure (corr.size() == 3);
        ensure (corr.rho(0, 0) == 1);
        ensure (corr.rho(0, 1) == 0.5);
        ensure (corr.rho(0, 2) == 0.4);
        ensure (corr.rho(1, 0) == 0.5);
        ensure (fabs(corr.rho(1, 1) - 1) <= eps);
        ensure (corr.rho(1, 2) == 0.5*0.4 + sqrt(1-0.5*0.5)*0.3);
        ensure (corr.rho(2, 0) == 0.4);
        ensure (corr.rho(2, 1) == 0.5*0.4 + sqrt(1-0.5*0.5)*0.3);
        ensure (fabs(corr.rho(2, 2) - 1) <= eps);
    }
}

template<class X>
void test_fms_swaption()
{
    std::vector<X> t(20), phi(20), sigma(20);
    X corr[] = {X(0.1),X(0.2),X(0.3)};
    auto freq = fms::fixed_income::frequency::semiannual;

    for (size_t i = 0; i < t.size(); ++i) {
        t[i] = (i + 1)/X(freq);
        phi[i] = X(0.05);
        sigma[i] = X(0.01);
    }

    fms::lmm<X,X> lmm(t.size(), t.data(), phi.data(), sigma.data(), fms::correlation<X>(t.size(), 3, corr));
    X pv = fms::swaption<X,X>(3, freq, X(0.05), 4, lmm, 1);
    pv = pv;
}

int main()
{
//    test_intB<double>();
    test_mean<double>();
    test_fms_correlation<double>();
    test_fms_brownian<double>();
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

    test_fms_pwflat<double>();
    //test_fms_pwflat<float>();

    test_fms_fixed_income_zero<double>();

    test_fms_pwflat_bootstrap<double>();

    test_fms_ho_lee<double>();

    test_fms_swaption<double>();
}