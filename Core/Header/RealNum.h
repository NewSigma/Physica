#ifndef _Physica_C_ReakNumber_H
#define _Physica_C_ReakNumber_H

#include "ComplexNum.h"

namespace Physica::Core {
    class RealNum : public ComplexNum {
    public:
        explicit RealNum(const Numerical& n);
        RealNum(const RealNum& n);
        explicit RealNum(const RealNum* n);
        ~RealNum() override;

        NumberType getType() const noexcept override;
        AbstractNum* operator+(const AbstractNum& n) const override;
        AbstractNum* operator-(const AbstractNum& n) const override;
        AbstractNum* operator*(const AbstractNum& n) const override;
        AbstractNum* operator/(const AbstractNum& n) const override;
        AbstractNum* operator^(const AbstractNum& n) const override;
        AbstractNum* operator-() const override;
        bool operator== (const AbstractNum& n) const override;
        bool operator> (const RealNum& n) const;
        bool operator< (const RealNum& n) const;
        bool operator>= (const RealNum& n) const;
        bool operator<= (const RealNum& n) const;
    };
}

#endif