// GR5260.cpp - test program
#include <cassert>
#include "root1d.h"
#include "fms_black.h"

using namespace fms;

void test_fms_black() {
    double f = 100, sigma = .2, k = 100, t = 0.25;

    double v = black::value(f, sigma, k, t);
    assert(fabs(v - 4.0) < .02);

    double s = sigma * sqrt(t);
    assert(v == black::value(f, s, k));
}

int main()
{
    test_fms_black();
}