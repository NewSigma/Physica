/*
 * Copyright 2022 WeiBo He.
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

namespace Physica::Core {
    template<class Derived>
    typename RValueVector<Derived>::RealType RValueVector<Derived>::norm() const {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    typename RValueVector<Derived>::RealType RValueVector<Derived>::squaredNorm() const {
        auto result = RealType::Zero();
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i).squaredNorm();
        return result;
    }

    template<class VectorType>
    typename RValueVector<VectorType>::ScalarType RValueVector<VectorType>::max() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result < temp)
                result = temp;
        }
        return result;
    }

    template<class VectorType>
    typename RValueVector<VectorType>::ScalarType RValueVector<VectorType>::min() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result > temp)
                result = temp;
        }
        return result;
    }

    template<class VectorType>
    typename RValueVector<VectorType>::ScalarType RValueVector<VectorType>::sum() const {
        assert(getLength() != 0);
        ScalarType result = ScalarType::Zero();
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i);
        return result;
    }

    template<class Derived>
    template<class OtherDerived>
    inline CrossProduct<Derived, OtherDerived>
    RValueVector<Derived>::crossProduct(const RValueVector<OtherDerived>& v) const noexcept {
        return CrossProduct(*this, v);
    }

    template<class Derived>
    template<class OtherDerived>
    typename RValueVector<Derived>::ScalarType
    RValueVector<Derived>::angleTo(const RValueVector<OtherDerived>& v) const noexcept {
        return arccos(Base::getDerived() * v.getDerived() / (norm() * v.norm()));
    }

    template<class VectorType>
    std::ostream& operator<<(std::ostream& os, const RValueVector<VectorType>& v) {
        os << '(';
        size_t length = v.getLength();
        if (length > 0) {
            --length;
            for (size_t i = 0; i < length; ++i)
                os << v.calc(i) << ", ";
            os << v.calc(length);
        }
        os << ')';
        return os;
    }
}
