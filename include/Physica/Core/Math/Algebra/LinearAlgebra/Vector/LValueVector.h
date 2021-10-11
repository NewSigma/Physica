/*
 * Copyright 2021 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "RValueVector.h"
#include "VectorBlock.h"

namespace Physica::Core {
    /**
     * \class LValueVector is base class of vectors that can be assigned to \class LValueVector
     * and other vectors can be assigned to this class.
     * In other words, you can take the address of elements in the vector.
     */
    template<class Derived>
    class LValueVector : public RValueVector<Derived> {
    public:
        using Base = RValueVector<Derived>;
        using typename Base::ScalarType;
    public:
        ~LValueVector() = default;
        /* Operators */
        LValueVector& operator=(const LValueVector& v);
        LValueVector& operator=(LValueVector&& v) { return (*this) = v; }
        template<class OtherVector>
        Derived& operator=(const RValueVector<OtherVector>& v);
        template<ScalarOption option, bool errorTrack>
        Derived& operator=(const Scalar<option, errorTrack>& s);
        [[nodiscard]] ScalarType& operator[](size_t index) { return Base::getDerived()[index]; }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return Base::getDerived()[index]; }
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueVector<OtherDerived>& v) const;
        [[nodiscard]] ScalarType calc(size_t index) const { return (*this)[index]; }
        VectorBlock<Derived> head(size_t to) { return VectorBlock<Derived>(Base::getDerived(), 0, to); }
        const VectorBlock<Derived> head(size_t to) const { return VectorBlock<Derived>(Base::getConstCastDerived(), 0, to); }
        VectorBlock<Derived> tail(size_t from) { return VectorBlock<Derived>(Base::getDerived(), from); }
        const VectorBlock<Derived> tail(size_t from) const { return VectorBlock<Derived>(Base::getConstCastDerived(), from); }
        VectorBlock<Derived> segment(size_t from, size_t to) { return VectorBlock<Derived>(Base::getDerived(), from, to); }
        const VectorBlock<Derived> segment(size_t from, size_t to) const { return VectorBlock<Derived>(Base::getConstCastDerived(), from, to); }
        /* Getters */
        [[nodiscard]] bool isZero() const;
    protected:
        LValueVector() = default;
        LValueVector(const LValueVector&) = default;
        LValueVector(LValueVector&&) noexcept = default;
    };

    template<class Derived>
    LValueVector<Derived>& LValueVector<Derived>::operator=(const LValueVector& v) {
        Base::getDerived().resize(v.getLength());
        v.assignTo(*this);
        return *this;
    }

    template<class Derived>
    template<class OtherVector>
    Derived& LValueVector<Derived>::operator=(const RValueVector<OtherVector>& v) {
        Base::getDerived().resize(v.getLength());
        v.assignTo(*this);
        return Base::getDerived();
    }

    template<class Derived>
    template<ScalarOption option, bool errorTrack>
    Derived& LValueVector<Derived>::operator=(const Scalar<option, errorTrack>& s) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = s;
        return Base::getDerived();
    }

    template<class Derived>
    template<class OtherDerived>
    void LValueVector<Derived>::assignTo(LValueVector<OtherDerived>& v) const {
        assert(v.getLength() == Base::getLength());
        for (size_t i = 0; i < Base::getLength(); ++i)
            v[i] = (*this)[i];
    }

    template<class Derived>
    bool LValueVector<Derived>::isZero() const {
        for(size_t i = 0; i < Base::getLength(); ++i)
            if (!(*this)[i].isZero())
                return false;
        return true;
    }

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const LValueVector<VectorType1>& v1, const LValueVector<VectorType2>& v2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type;
        const auto len = v1.getLength();
        assert(len == v2.getLength());
        auto result = ResultType::Zero();
        for(size_t i = 0; i < len; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    template<class Derived, class OtherDerived>
    void operator+=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        for (size_t i = 0; i < v1.getLength(); ++i)
            v1[i] = v1[i] + v2.calc(i);
    }

    template<class Derived, class OtherDerived>
    void operator-=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        for (size_t i = 0; i < v1.getLength(); ++i)
            v1[i] = v1[i] - v2.calc(i);
    }

    template<class VectorType, class ScalarType>
    inline void operator+=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v + n; }

    template<class VectorType, class ScalarType>
    inline void operator-=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v - n; }

    template<class VectorType, class ScalarType>
    inline void operator*=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v * n; }

    template<class VectorType, class ScalarType>
    inline void operator/=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v / n; }

    template<class VectorType>
    std::ostream& operator<<(std::ostream& os, const LValueVector<VectorType>& v) {
        os << '(';
        size_t length = v.getLength();
        if (length) {
            --length;
            for (size_t i = 0; i < length; ++i)
                os << v[i] << ", ";
            os << v[length];
        }
        os << ')';
        return os;
    }
}
