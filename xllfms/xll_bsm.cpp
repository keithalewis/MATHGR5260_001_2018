// xll_bsm.cpp - Black-Scholes/Merton model add-in.
#include "../GR5260/fms_bsm.h"
#include "../xll12/xll/xll.h"

using namespace xll;

AddIn xai_bsm_value(
    Function(XLL_DOUBLE, L"?xll_bsm_value", L"BSM.VALUE")
    .Arg(XLL_DOUBLE, L"r", L"is the risk-neutral continuously compounded interest rate.")
    .Arg(XLL_DOUBLE, L"s", L"is the spot rate.")
    .Arg(XLL_DOUBLE, L"sigma", L"is the annual volatility.")
    .Arg(XLL_DOUBLE, L"k", L"is the strike.")
    .Arg(XLL_DOUBLE, L"t", L"is the time in years to expiration.")
    .Category(L"GR5260")
    .FunctionHelp(L"Return the Black put option forward value.")
    .Documentation(L"TBA")
);
double WINAPI xll_bsm_value(double r, double s, double sigma, double k, double t)
{
#pragma XLLEXPORT
    double result = std::numeric_limits<double>::quiet_NaN();

    try {
        result = fms::bsm::value(r, s, sigma, k, t);
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return result;
}

AddIn xai_bsm_delta(
    Function(XLL_DOUBLE, L"?xll_bsm_delta", L"BSM.DELTA")
    .Arg(XLL_DOUBLE, L"r", L"is the risk-neutral continuously compounded interest rate.")
    .Arg(XLL_DOUBLE, L"s", L"is the spot rate.")
    .Arg(XLL_DOUBLE, L"sigma", L"is the annual volatility.")
    .Arg(XLL_DOUBLE, L"k", L"is the strike.")
    .Arg(XLL_DOUBLE, L"t", L"is the time in years to expiration.")
    .Category(L"GR5260")
    .FunctionHelp(L"Return the Black put option forward delta.")
    .Documentation(L"TBA")
);
double WINAPI xll_bsm_delta(double r, double s, double sigma, double k, double t)
{
#pragma XLLEXPORT
    double result = std::numeric_limits<double>::quiet_NaN();

    try {
        result = fms::bsm::delta(r, s, sigma, k, t);
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return result;
}