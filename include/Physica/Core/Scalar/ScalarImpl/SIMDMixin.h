/*
 * Copyright 2024-2026 Weibo He.
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

#include <memory>
#include <utility>
#include "Physica/CRTP.h"
#include "Physica/Core/Math/Algebra/Canonicalization.h"

namespace Physica {
    template<class Derived>
    class SIMDMixin : public CRTP<SIMDMixin<Derived>> {
        using This = SIMDMixin<Derived>;
        using Base = CRTP<This>;
        using TraitsType = Traits<Derived>;
    public:
        constexpr static bool isSeparatable = TraitsType::isSeparatable;

        using ScalarType = TraitsType::ScalarType;
        using ValueType = TraitsType::ValueType;
        using GradType = TraitsType::GradType;
        using FullRealType = TraitsType::FullRealType;
        using RealType = TraitsType::RealType;
        using BoolSIMDType = TraitsType::BoolSIMDType;
        using MachineType = TraitsType::MachineType;

        class Iterator;
    public:
        constexpr ~SIMDMixin() = default;
        /* Operations */
        [[nodiscard]] FullRealType squaredNorm() const;
        [[nodiscard]] FullRealType swapRealImag() const;
        [[nodiscard]] FullRealType gatherRealImag() const noexcept;
        [[nodiscard]] FullRealType scatterRealImag() const noexcept;
        /* Getters */
        [[nodiscard, gnu::always_inline]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto end(this auto&&) noexcept;
        [[nodiscard]] constexpr static int size() noexcept { return TraitsType::Size; }
        [[nodiscard]] ValueType value() const noexcept;
        [[nodiscard]] FullRealType asReal() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept { return ScalarType::isComplex(); }
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept { return ScalarType::isDiffable(); }
    protected:
        constexpr SIMDMixin() = default;
        constexpr SIMDMixin(const This&) = default;
        constexpr SIMDMixin(This&&) = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
    };

    template<class Derived>
    auto SIMDMixin<Derived>::squaredNorm() const -> FullRealType {
        const FullRealType x2 = square(asReal());
        if constexpr (isComplex())
            return x2 + x2.swapRealImag();
        else
            return x2;
    }

    template<class Derived>
    auto SIMDMixin<Derived>::swapRealImag() const -> FullRealType {
        const auto x = asReal();
        if constexpr (ScalarType::Prec == Float32)
            return x.template shuffle<1, 0, 3, 2>();
        else {
            constexpr int Size1 = isComplex() ? size() * 2 : size();
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
    auto SIMDMixin<Derived>::gatherRealImag() const noexcept -> FullRealType {
        const auto x = asReal();
        constexpr int Size1 = isComplex() ? size() * 2 : size();
        if constexpr (Size1 == 2)
            return x;
        if constexpr (Size1 == 4)
            return x.template permute<0, 2, 1, 3>();
        else if constexpr (Size1 == 8)
            return x.template permute<0, 2, 4, 6, 1, 3, 5, 7>();
        else {
            static_assert(Size1 == 16, "[Error]: Unexpected size");
            return x.template permute<0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15>();
        }
    }

    template<class Derived>
    auto SIMDMixin<Derived>::scatterRealImag() const noexcept -> FullRealType {
        const auto x = asReal();
        constexpr int Size1 = isComplex() ? size() * 2 : size();
        if constexpr (Size1 == 2)
            return x;
        if constexpr (Size1 == 4)
            return x.template permute<0, 2, 1, 3>();
        else if constexpr (Size1 == 8)
            return x.template permute<0, 4, 1, 5, 2, 6, 3, 7>();
        else {
            static_assert(Size1 == 16, "[Error]: Unexpected size");
            return x.template permute<0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15>();
        }
    }

    template<class Derived>
    auto SIMDMixin<Derived>::value() const noexcept -> ValueType {
        if constexpr (isDiffable())
            return Base::getDerived_host().value();
        else
            return Base::getDerived_host();
    }

    template<class Derived>
    auto SIMDMixin<Derived>::asReal() const noexcept -> FullRealType {
        if constexpr (isComplex())
            return Base::getDerived_host().asReal();
        else
            return Base::getDerived_host();
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::begin(this auto&& self) noexcept {
        return Iterator(std::addressof(self), 0);
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::end(this auto&& self) noexcept {
        return Iterator(std::addressof(self), TraitsType::Size);
    }

    template<class Derived>
    class SIMDMixin<Derived>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::random_access_iterator_tag;
        using difference_type = int;
        using value_type = typename SIMDMixin::ScalarType;
        using reference = const value_type;
        using const_reference = const value_type;
    private:
        const Derived* pack = nullptr;
        difference_type index{};
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] constexpr Iterator(const Derived* pack, difference_type index) noexcept;
        constexpr Iterator(const This&) = default;
        constexpr Iterator(This&&) noexcept = default;
        constexpr ~Iterator() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        [[gnu::always_inline]] constexpr This& operator++() noexcept;
        [[gnu::always_inline]] constexpr This& operator--() noexcept;
        [[gnu::always_inline]] constexpr This& operator+=(difference_type n) noexcept;
        [[gnu::always_inline]] constexpr This& operator-=(difference_type n) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator++(int) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator--(int) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr reference operator*() const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr reference operator[](difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto operator<=>(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator+(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator-(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr difference_type operator-(const This& other) const noexcept;
        /* Friends */
        [[gnu::always_inline]] friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<class Derived>
    constexpr SIMDMixin<Derived>::Iterator::Iterator(const Derived* pack, difference_type index) noexcept : pack(pack), index(index) {}

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator++() noexcept -> This& {
        index += 1;
        return *this;
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator--() noexcept -> This& {
        index -= 1;
        return *this;
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator+=(difference_type n) noexcept -> This& {
        index += n;
        return *this;
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator-=(difference_type n) noexcept -> This& {
        index -= n;
        return *this;
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(pack, index + 1));
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, This(pack, index - 1));
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator*() const noexcept -> reference {
        return (*pack)[index];
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator[](difference_type n) const noexcept -> reference {
        return (*pack)[index + n];
    }

    template<class Derived>
    constexpr bool SIMDMixin<Derived>::Iterator::operator==(const This& other) const noexcept {
        return index == other.index;
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator<=>(const This& other) const noexcept {
        return index <=> other.index;
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(pack, index + n);
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(pack, index - n);
    }

    template<class Derived>
    constexpr auto SIMDMixin<Derived>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return index - other.index;
    }

    [[nodiscard]] auto operator*(const Scalar auto& x, const Packet auto p) noexcept {
        return p * x;
    }

    [[nodiscard]] auto operator*(const Packet auto p1, const Packet auto p2) noexcept {
        if constexpr (canonicalized(p1, p2))
            return p1.operator*(p2);
        else
            return p2 * p1;
    }
}

namespace Physica {
    template<class T>
    class Traits<SIMDMixin<T>> {
    public:
        using Derived = T;
    };
}
