// root1d.h - One dimensional root finding.
#pragma once
#include <functional>

namespace fms::root1d {

    // https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-Virtual_Interface
    template<class X = double>
    class abstract_solver {
    public:
        // next guess at the root
        X next()
        {
            return _next();
        }
        // no more iterations needed
        bool done()
        {
            return _done();
        }
        X solve()
        {
            X x;

            do {
                x = next();
            } while (!done());

            return x;
        }
    private:
        virtual X _next() = 0;
        virtual bool _done() = 0;
    };

} // fms::root1d
