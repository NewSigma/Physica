/*
 * Copyright 2020-2024 Weibo He.
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

#include <cstddef>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Vector.h"

namespace Physica::Core {
    template<Scalar T, size_t Length = Dynamic>
    class UnitVector : public RValueVector<UnitVector<T, Length>> {
        using This = UnitVector<T, Length>;
        using Base = RValueVector<This>;

        size_t index1;
        size_t length;
    public:
        UnitVector(size_t index1_, size_t length_);
        UnitVector(const This&) = default;
        UnitVector(This&&) noexcept = default;
        ~UnitVector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return length; }
    };

    template<Scalar T, size_t Length>
    UnitVector<T, Length>::UnitVector(size_t index1_, size_t length_) : index1(index1_), length(length_) {}

    template<Scalar T, size_t Length>
    void UnitVector<T, Length>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(index1, obj.index1);
        std::swap(length, obj.length);
    }

    template<Scalar T, size_t Length>
    T UnitVector<T, Length>::calc(size_t index) const {
        return index == index1 ? T(1) : T(0);
    }
}

namespace Physica {
    template<Scalar T, size_t Length>
    class Traits<UnitVector<T, Length>> {
    public:
        using ScalarType = T;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
