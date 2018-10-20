// fms_fixed_income_instrument.h - Fixed Income instruments interface class
#pragma once

namespace fms::fixed_income {

    enum class frequency {
        annual = 1,
        semiannual = 2,
        quarterly = 4,
        monthly = 12,
    };

    // NVI instrument interface.
    template<class U = double, class C = double>
    class instrument {
        size_t n;
        U* u;
        C* c;
    public:
        // For pre-allocated and unmanaged memory.
        instrument(size_t n = 0, U* u = 0, C* c = 0)
            : n(n), u(u), c(c)
        { }
        virtual ~instrument() {}
        // convenience functions
        bool operator==(const instrument& i) const
        {
            return size() == i.size()
                && std::equal(time(), time() + size(), i.time())
                && std::equal(cash(), cash() + size(), i.cash());
        }
        bool operator!=(const instrument& i) const
        {
            return !operator==(i);
        }

        size_t   size() const { return _size(); }
        const U* time() const { return _time(); }
        const C* cash() const { return _cash(); }

        // time of last cash flow
        U termination() const 
        {
            if (size() == 0)
                return std::numeric_limits<U>::quite_NaN();

            return time()[size()-1];
        }
    private:
        // Override in derived classes.
        virtual size_t _size() const
        {
            return n;
        }
        virtual const U* _time() const
        {
            return u;
        }
        virtual const C* _cash() const
        {
            return c;
        }
    };

} // fms::fixed_income
