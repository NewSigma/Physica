/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POLYNOMIAL_H
#define PHYSICA_POLYNOMIAL_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Utils/CStyleArray/CStyleArray.h"

namespace Physica::Core {
    /*!
     * Class Polynomial provides operations around polynomials.
     *
     * A polynomial is:
     * y(x) = a0 + a1 x + a2 x ^ 2 + ... + an x ^ n
     */
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class Polynomial {
    public:
        //data stores the coeficients of the polynomial.
        CStyleArray<Scalar<type, errorTrack>, Dynamic> data;
    public:
        Polynomial() = default;
        Polynomial(const Polynomial& p) = default;
        Polynomial(Polynomial&& p) noexcept : data(std::move(p.data)) {}
        ~Polynomial() = default;
        /* Operators */
        Polynomial& operator=(const Polynomial& p) = default;
        Polynomial& operator=(Polynomial&& p) noexcept { data = std::move(p.data); return *this; }
        Scalar<type, errorTrack> operator()(const Scalar<type, errorTrack>& s) const;
    };

    /*!
     * Returns the value of this polynomial when x = s.
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> Polynomial<type, errorTrack>::operator()(const Scalar<type, errorTrack>& s) const {
        if(data.empty())
            return Scalar<type, errorTrack>::getZero();
        Scalar<type, errorTrack> result(data[0]);
        Scalar<type, errorTrack> temp(s);
        const auto length = data.getLength();
        for(size_t i = 1; i < length; ++i) {
            result += temp * data[i];
            temp *= s;
        }
        return result;
    }
}

#endif
