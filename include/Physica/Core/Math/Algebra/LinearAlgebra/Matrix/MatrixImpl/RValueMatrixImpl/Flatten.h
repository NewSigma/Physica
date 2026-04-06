/*
 * Copyright 2022-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class FlattenR<M> : public RValueVector<FlattenR<M>> {
        using This = FlattenR<M>;
        using Base = RValueVector<FlattenR<M>>;

        LazyDestroy<M> mat;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        FlattenR(M&& mat_) : mat(std::forward<M>(mat_)) {}
        FlattenR(const This&) = default;
        FlattenR(This&&) noexcept = default;
        ~FlattenR() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        [[nodiscard]] auto grads(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
    };

    template<Matrix M>
    auto FlattenR<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().flatten();
    }

    template<Matrix M>
    auto FlattenR<M>::grads(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).grads().flatten();
    }

    template<Matrix M>
    auto FlattenR<M>::calc(size_t index) const -> T {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor);
    }

    template<Matrix M>
    auto FlattenR<M>::calc_value(size_t index) const -> Tv {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calc_value(mat.rowFromMajorMinor(major, minor),
                              mat.colFromMajorMinor(major, minor));
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<FlattenR<M>> {
        using M1 = std::remove_reference_t<M>;
    public:
        using ScalarType = M1::ScalarType;
        constexpr static size_t SizeAtCompile = M1::RowAtCompile * M1::ColAtCompile;

        constexpr static bool FastAssign = false;
    };
}
