/*
 * Copyright 2023-2026 Weibo He.
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

namespace Physica {
    template<Scalar T, int Size>
    class BoolSIMD : private Traits<BoolSIMD<T, Size>>::Pack {
        using This = BoolSIMD<T, Size>;
        using Base = Traits<This>::Pack;
        using HalfType = std::conditional<sizeof(Base) * CHAR_BIT != 128, BoolSIMD<T, Size / 2>, Empty>::type;
    public:
        BoolSIMD() = default;
        explicit BoolSIMD(Base value) : Base(value) {}
        using Base::Base;
        BoolSIMD(const This&) = default;
        BoolSIMD(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] This operator&&(This other) const noexcept;
        [[nodiscard]] This operator||(This other) const noexcept;
        [[nodiscard]] This operator^(This other) const noexcept;
        [[nodiscard]] This operator!() const noexcept;
        /* Operations */
        [[nodiscard]] bool horizontal_and() const;
        [[nodiscard]] bool horizontal_or() const;
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] constexpr static size_t getSize() { return Size; }
        [[nodiscard]] Base& toMachine() noexcept { return *this; }
        [[nodiscard]] const Base& toMachine() const noexcept { return *this; }
        [[nodiscard]] auto getLow() const noexcept { return HalfType(Base::get_low()); }
        [[nodiscard]] auto getHigh() const noexcept { return HalfType(Base::get_high()); }
    };

    template<Scalar T, int Size>
    auto BoolSIMD<T, Size>::operator&&(This other) const noexcept -> This {
        return This(toMachine() && other.toMachine());
    }

    template<Scalar T, int Size>
    auto BoolSIMD<T, Size>::operator||(This other) const noexcept -> This {
        return This(toMachine() && other.toMachine());
    }

    template<Scalar T, int Size>
    auto BoolSIMD<T, Size>::operator^(This other) const noexcept -> This {
        return This(toMachine() ^ other.toMachine());
    }

    template<Scalar T, int Size>
    auto BoolSIMD<T, Size>::operator!() const noexcept -> This {
        return This(!toMachine());
    }

    template<Scalar T, int Size>
    bool BoolSIMD<T, Size>::horizontal_and() const {
        return Physica::horizontal_and(toMachine());
    }

    template<Scalar T, int Size>
    bool BoolSIMD<T, Size>::horizontal_or() const {
        return Physica::horizontal_or(toMachine());
    }

    template<Scalar T>
    class BoolSIMD<T, 1> {
        using This = BoolSIMD<T, 1>;
        using Base = Traits<This>::Pack;
    private:
        bool b;
    public:
        BoolSIMD() = default;
        explicit BoolSIMD(bool value) : b(value) {}
        BoolSIMD(const This&) = default;
        BoolSIMD(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] bool horizontal_or() const { return b; }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return 1; }
        [[nodiscard]] constexpr static size_t getSize() { return 1; }
    };
}

namespace Physica {
    template<Scalar T, int Size>
    class Traits<BoolSIMD<T, Size>> {
    public:
        using ScalarType = T;
    private:
        constexpr static bool isSinglePrec = ScalarType::Prec == Float;

        using Size2Type = std::conditional<isSinglePrec, void, Vec2db>::type;
        using Size4Type = std::conditional<isSinglePrec, Vec4fb, Vec4db>::type;
        using Size8Type = std::conditional<isSinglePrec, Vec8fb, Vec8db>::type;
        using Size16Type = std::conditional<isSinglePrec, Vec16fb, void>::type;
        using Type1 = std::conditional<Size == 2, Size2Type, Size4Type>::type;
        using Type2 = std::conditional<Size == 8, Size8Type, Size16Type>::type;
    public:
        using Pack = std::conditional<Size <= 4, Type1, Type2>::type;
    };
}
