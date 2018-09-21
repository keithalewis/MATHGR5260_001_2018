// xll_black.cpp - Black model add-in.
#include "../GR5260/fms_black.h"
#include "../xll12/xll/xll.h"

using namespace xll;

AddIn xai_black_value(
    Function(XLL_DOUBLE, L"?xll_black_value", L"BLACK.VALUE")
    .Arg(XLL_DOUBLE, L"f", L"is ???")
    .Arg(XLL_DOUBLE, L"sigma", L"is ???")
    .Arg(XLL_DOUBLE, L"k", L"is ???")
    .Arg(XLL_DOUBLE, L"t", L"is ???")
    .Category(L"GR5260")
    .FunctionHelp(L"Return the Black put option forward value.")
    .Documentation(L"TBA")
);
double WINAPI xll_black_value(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
    double result = std::numeric_limits<double>::quiet_NaN();

    try {
        result = fms::black::value(f, sigma, k, t);
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return result;
}

AddIn xai_black_delta(
    Function(XLL_DOUBLE, L"?xll_black_delta", L"BLACK.DELTA")
    .Arg(XLL_DOUBLE, L"f", L"is ???")
    .Arg(XLL_DOUBLE, L"sigma", L"is ???")
    .Arg(XLL_DOUBLE, L"k", L"is ???")
    .Arg(XLL_DOUBLE, L"t", L"is ???")
    .Category(L"GR5260")
    .FunctionHelp(L"Return the Black put option forward delta.")
    .Documentation(L"TBA")
);
double WINAPI xll_black_delta(double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
    double result = std::numeric_limits<double>::quiet_NaN();

    try {
        result = fms::black::delta(f, sigma, k, t);
    }
    catch(const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return result;
}