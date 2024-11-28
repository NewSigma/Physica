/*
 * Copyright 2024 Weibo He.
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

#include "Physica/CRTPBase.h"
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<class Derived>
    class SIMDBase : public CRTPBase<SIMDBase<Derived>> {
        using This = SIMDBase<Derived>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;
    public:
        constexpr static int Size = TraitsType::Size;

        using ScalarType = TraitsType::ScalarType;
        using RealType = TraitsType::RealType;
        using FullRealType = TraitsType::FullRealType;
        
        constexpr static bool isComplex = ScalarType::isComplex;
        static_assert(isComplex != std::is_same<RealType, FullRealType>::value, "[Error]: Inconsistent declaration");
    public:
        constexpr ~SIMDBase() = default;
        /* Operations */
        [[nodiscard]] FullRealType squaredNorm() const;
        [[nodiscard]] FullRealType swapRealImag() const;
        /* Getters */
        [[nodiscard]] FullRealType asReal() const;
    protected:
        constexpr SIMDBase() = default;
        constexpr SIMDBase(const This&) = default;
        constexpr SIMDBase(This&&) = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
    };

    template<class Derived>
    SIMDBase<Derived>::FullRealType SIMDBase<Derived>::squaredNorm() const {
        const FullRealType x2 = square(asReal());
        if constexpr (isComplex)
            return x2 + x2.swapRealImag();
        else
            return x2;
    }

    template<class Derived>
    SIMDBase<Derived>::FullRealType SIMDBase<Derived>::swapRealImag() const {
        const auto x = asReal();
        if constexpr (ScalarType::Option == Float32)
            return x.template shuffle<1, 0, 3, 2>();
        else {
            constexpr int Size1 = isComplex ? Size * 2 : Size;
            if constexpr (Size1 == 2)
                return x.template shuffle<1, 0>();
            else if constexpr (Size1 == 4)
                return x.template shuffle<1, 0, 1, 0>();
            else {
                static_assert(Size1 == 8, "[Error]: Unexpected size");
                return x.template shuffle<1, 0, 1, 0, 1, 0, 1, 0>();
            }
        }
    }

    template<class Derived>
    SIMDBase<Derived>::FullRealType SIMDBase<Derived>::asReal() const {
        if constexpr (isComplex)
            return Base::getDerived().asReal();
        else
            return Base::getDerived();
    }
}

namespace Physica {
    template<class T>
    class Traits<SIMDBase<T>> {
    public:
        using Derived = T;
    };
}
