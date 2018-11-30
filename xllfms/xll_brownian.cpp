// xll_brownian.cpp
#include <random>
#include "../GR5260/fms_brownian.h"
#include "xllfms.h"

using namespace fms;
using namespace xll;

AddIn xai_fms_brownian(
    Function(XLL_HANDLE, L"?xll_fms_brownian", L"FMS.BROWNIAN")
    .Arg(XLL_HANDLE, L"correlation", L"is a handle to a fms::correlation object")
    .Uncalced()
    .Category(L"GR5260")
    .FunctionHelp(L"Return a handle to a fms::brownian object.")
    .Documentation(L"TBD")
);
HANDLEX WINAPI xll_fms_brownian(HANDLEX c)
{
#pragma XLLEXPORT
    handlex h;

    try {
        handle<fms::correlation<>> c_(c);
        ensure (c_);
        handle<fms::brownian<>> h_(new fms::brownian<>(*c_));
        ensure (h_);
        h = h_.get();
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return h;
}

AddIn xai_fms_brownian_reset(
    Function(XLL_HANDLE, L"?xll_fms_brownian_reset", L"FMS.BROWNIAN.RESET")
    .Arg(XLL_HANDLE, L"handle", L"is a handle to a fms::brownain object")
    .Category(L"GR5260")
    .FunctionHelp(L"Reset samples to time 0.")
    .Documentation(L"TBD")
);
HANDLEX WINAPI xll_fms_brownian_reset(HANDLEX b)
{
#pragma XLLEXPORT
    try {
        handle<fms::brownian<>> b_(b);
        ensure (b_);
        b_->reset();
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return b;
}

AddIn xai_fms_brownian_advance(
    Function(XLL_HANDLE, L"?xll_fms_brownian_advance", L"FMS.BROWNIAN.ADVANCE")
    .Arg(XLL_HANDLE, L"handle", L"is a handle to a fms::brownian object")
    .Arg(XLL_DOUBLE, L"time", L"is the time to advance the brownian object.")
    .Category(L"GR5260")
    .FunctionHelp(L"Advance a Brownian sample to a give time.")
    .Documentation(L"TBD")
);
HANDLEX WINAPI xll_fms_brownian_advance(HANDLEX b, double t)
{
#pragma XLLEXPORT
    try {
        handle<fms::brownian<>> b_(b);
        ensure (b_);
        b_->advance(t, fms::dre);
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return b;
}

AddIn xai_fms_brownian_sample(
    Function(XLL_FP, L"?xll_fms_brownian_sample", L"FMS.BROWNIAN.SAMPLE")
    .Arg(XLL_HANDLE, L"handle", L"is a handle to a fms::brownain object")
    .Arg(XLL_FP, L"times", L"is an array of increasing times.")
    .Volatile()
    .Category(L"GR5260")
    .FunctionHelp(L"Return samples of correlated Brownian motion.")
    .Documentation(L"TBD")
);
_FP12* WINAPI xll_fms_brownian_sample(HANDLEX b, _FP12* pt)
{
#pragma XLLEXPORT
    static xll::FP12 result;

    try {
        handle<fms::brownian<>> b_(b);
        ensure (b_);
        result.resize(size(*pt), b_->size());
        b_->reset();
        for (INT32 i = 0; i < size(*pt); ++i) {
            b_->advance(pt->array[i], fms::dre);
            std::copy(b_->data(), b_->data() + b_->size(),
                result.array() + i * result.columns());
        }
    }
    catch (const std::exception& ex) {
        XLL_ERROR(ex.what());
    }

    return result.get();
}
