/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    template<class Derived>
    class SIMDBase : public CRTPBase<SIMDBase<Derived>> {
        using This = SIMDBase<Derived>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;
    public:
        constexpr static int Size = TraitsType::Size;
        constexpr static bool isSeparatable = TraitsType::isSeparatable;

        using ScalarType = TraitsType::ScalarType;
        using ValueType = TraitsType::ValueType;
        using GradType = TraitsType::GradType;
        using FullRealType = TraitsType::FullRealType;
        using RealType = TraitsType::RealType;
        using BoolSIMDType = TraitsType::BoolSIMDType;
        using MachineType = TraitsType::MachineType;
        
        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static bool isDiffable = ScalarType::isDiffable;
    public:
        constexpr ~SIMDBase() = default;
        /* Operations */
        [[nodiscard]] inline FullRealType squaredNorm() const;
        [[nodiscard]] inline FullRealType swapRealImag() const;
        [[nodiscard]] inline FullRealType permRealImag() const noexcept;
        [[nodiscard]] inline FullRealType scatterRealImag() const noexcept;
        /* Getters */
        [[nodiscard]] constexpr static int size() { return Size; }
        [[nodiscard]] inline ValueType value() const;
        [[nodiscard]] inline FullRealType asReal() const;
    protected:
        constexpr SIMDBase() = default;
        constexpr SIMDBase(const This&) = default;
        constexpr SIMDBase(This&&) = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
    };

    template<class Derived>
    inline auto SIMDBase<Derived>::squaredNorm() const -> FullRealType {
        const FullRealType x2 = square(asReal());
        if constexpr (isComplex)
            return x2 + x2.swapRealImag();
        else
            return x2;
    }

    template<class Derived>
    inline auto SIMDBase<Derived>::swapRealImag() const -> FullRealType {
        const auto x = asReal();
        if constexpr (ScalarType::Prec == Float32)
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
    inline auto SIMDBase<Derived>::permRealImag() const noexcept -> FullRealType {
        const auto x = asReal();
        constexpr int Size1 = isComplex ? Size * 2 : Size;
        if constexpr (Size1 == 2)
            return x;
        if constexpr (Size1 == 4)
            return x.template permute<0, 2, 1, 3>();
        else if constexpr (Size1 == 8)
            return x.template permute<0, 2, 4, 6, 1, 3, 5, 7>();
        else {
            static_assert(Size == 16, "[Error]: Unexpected size");
            return x.template permute<0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15>();
        }
    }

    template<class Derived>
    inline auto SIMDBase<Derived>::scatterRealImag() const noexcept -> FullRealType {
        const auto x = asReal();
        constexpr int Size1 = isComplex ? Size * 2 : Size;
        if constexpr (Size1 == 2)
            return x;
        if constexpr (Size1 == 4)
            return x.template permute<0, 2, 1, 3>();
        else if constexpr (Size1 == 8)
            return x.template permute<0, 4, 1, 5, 2, 6, 3, 7>();
        else {
            static_assert(Size == 16, "[Error]: Unexpected size");
            return x.template permute<0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15>();
        }
    }

    template<class Derived>
    inline auto SIMDBase<Derived>::value() const -> ValueType {
        if constexpr (isDiffable)
            return Base::getDerived_host().value();
        else
            return Base::getDerived_host();
    }

    template<class Derived>
    inline auto SIMDBase<Derived>::asReal() const -> FullRealType {
        if constexpr (isComplex)
            return Base::getDerived_host().asReal();
        else
            return Base::getDerived_host();
    }
}

namespace Physica {
    template<class T>
    class Traits<SIMDBase<T>> {
    public:
        using Derived = T;
    };
}
