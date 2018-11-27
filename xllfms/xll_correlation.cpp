// xll_correlation.cpp
#include "../GR5260/fms_correlation.h"
#include "../xll12/xll/xll.h"

using namespace fms;
using namespace xll;

//??? Implement the add-in FMS.CORRELATION that takes an n x d matrix in lower triangular form.
//??? It should return a xll::handle<fms::correlation>.
AddIn xai_fms_correlation(
    Function(XLL_HANDLE, L"?xll_fms_correlation", L"FMS.CORRELATION")
    .Arg(XLL_FP, L"matrix", L"is an n x d matix in lower triangular form.")
    .Uncalced()
    .Category(L"GR5260")
    .FunctionHelp(L"Return handle to correlation object.")
    .Documentation(L"TBA")
);
HANDLEX WINAPI xll_fms_correlation(_FP12* pm)
{
#pragma XLLEXPORT
    handlex h;

    try {

        using layout = fms::correlation<>::layout;
        xll::handle<fms::correlation<>> h_(new fms::correlation<>(pm->rows + 1, pm->columns + 1, pm->array, layout::lower));
        h = h_.get();
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return h;
}

//??? Implment FMS.CORRELATION.RHO that calls fms::correlation::rho.
AddIn xai_fms_correlation_rho(
    Function(XLL_DOUBLE, L"?xll_fms_correlation_rho", L"FMS.CORRELATION.RHO")
    .Arg(XLL_HANDLE, L"handle", L"is handle to an fms::correlation object.")
    .Arg(XLL_WORD, L"i", L"is the row index.")
    .Arg(XLL_WORD, L"j", L"is the column index.")
    .Category(L"GR5260")
    .FunctionHelp(L"Return i,j correlation.")
    .Documentation(L"TBA")
);
HANDLEX WINAPI xll_fms_correlation_rho(HANDLEX h, WORD i, WORD j)
{
#pragma XLLEXPORT
    double rho = std::numeric_limits<double>::quiet_NaN();

    try {
        // i , j < n
        xll::handle<fms::correlation<>> h_(h);
        // ensure (h_);
        rho = h_->rho(i, j);
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return rho;
}
