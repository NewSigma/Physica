/*
 * Copyright 2023-2024 Weibo He.
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
    template<Scalar T, size_t Size>
    class BoolSIMD : private Traits<BoolSIMD<T, Size>>::BaseType {
        using This = BoolSIMD<T, Size>;
    public:
        using Base = typename Traits<This>::BaseType;
    public:
        BoolSIMD() = default;
        explicit BoolSIMD(Base value) : Base(value) {}
        using Base::Base;
        BoolSIMD(const BoolSIMD&) = default;
        BoolSIMD(BoolSIMD&&) noexcept = default;
        /* Operators */
        BoolSIMD& operator=(const BoolSIMD&) = default;
        BoolSIMD& operator=(BoolSIMD&&) noexcept = default;
        /* Operations */
        [[nodiscard]] inline bool horizontal_or() const;
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] constexpr static size_t getSize() { return Size; }
        [[nodiscard]] Base& toMachine() noexcept { return *this; }
        [[nodiscard]] const Base& toMachine() const noexcept { return *this; }
    };

    template<Scalar T, size_t Size>
    [[nodiscard]] inline bool BoolSIMD<T, Size>::horizontal_or() const {
        return Physica::horizontal_or(toMachine());
    }

    template<Scalar T>
    class BoolSIMD<T, 1> {
        using This = BoolSIMD<T, 1>;
    public:
        using Base = typename Traits<This>::BaseType;
    private:
        bool b;
    public:
        BoolSIMD() = default;
        explicit BoolSIMD(bool value) : b(value) {}
        BoolSIMD(const BoolSIMD&) = default;
        BoolSIMD(BoolSIMD&&) noexcept = default;
        /* Operators */
        BoolSIMD& operator=(const BoolSIMD&) = default;
        BoolSIMD& operator=(BoolSIMD&&) noexcept = default;
        /* Operations */
        [[nodiscard]] inline bool horizontal_or() const { return b; }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return 1; }
        [[nodiscard]] constexpr static size_t getSize() { return 1; }
    };
}

namespace Physica {
    template<Scalar T, size_t Size>
    class Traits<BoolSIMD<T, Size>> {
    public:
        using ScalarType = T;
    private:
        constexpr static bool isSinglePrec = ScalarType::Option == Float;
            
        using Size2Type = typename std::conditional<isSinglePrec, void, Vec2db>::type;
        using Size4Type = typename std::conditional<isSinglePrec, Vec4fb, Vec4db>::type;
        using Size8Type = typename std::conditional<isSinglePrec, Vec8fb, Vec8db>::type;
        using Size16Type = typename std::conditional<isSinglePrec, Vec16fb, void>::type;
        using Type1 = typename std::conditional<Size == 2, Size2Type, Size4Type>::type;
        using Type2 = typename std::conditional<Size == 8, Size8Type, Size16Type>::type;
    public:
        using BaseType = typename std::conditional<Size <= 4, Type1, Type2>::type;
    };
}
