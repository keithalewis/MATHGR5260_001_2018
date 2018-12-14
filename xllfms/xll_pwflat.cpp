// xll_pwflat.cpp - Piecewise flat curves.
#include "../GR5260/fms_pwflat.h"
#include "../xll12/xll/xll.h"

using namespace xll;

class AddIn xai_fms_pwflat_value(
    Function(XLL_DOUBLE, L"?xll_fms_pwflat_value", L"PWFLAT.VALUE")
    .Arg(XLL_DOUBLE, L"u", L"is the time at which to value the curve.")
    .Arg(XLL_FP, L"t", L"is a list of times.")
    .Arg(XLL_FP, L"f", L"is a list of forwards.")
    .Arg(XLL_DOUBLE, L"_f", L"is an optional forward to extrapolate the curve.")
    .Category(L"GR5260")
    .FunctionHelp(L"Return the value of a piecewise flat curve.")
    .Documentation(fms_pwflat_doc)
);
double WINAPI xll_fms_pwflat_value(double u, _FP12* pt, _FP12* pf, double _f)
{
#pragma XLLEXPORT
    double result;

    try {
        ensure (size(*pt) == size(*pf));
        result = fms::pwflat::value(u, size(*pt), pt->array, pf->array, _f);
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());

        result = std::numeric_limits<double>::quiet_NaN();
    }

    return result;
}
