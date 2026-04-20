/*
 * Copyright 2026 Weibo He.
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

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class MinorDiagR final : public RValueVector<MinorDiagR<M>> {
        using This = MinorDiagR<M>;
        using Base = RValueVector<This>;
    private:
        LazyDestroy<M> mat;
        ssize_t shift;
    public:
        MinorDiagR(M mat, ssize_t shift);
        MinorDiagR(const This&) = default;
        MinorDiagR(This&&) = default;
        ~MinorDiagR() = default;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(size_t index) const noexcept;
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() - std::abs(shift); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M>
    MinorDiagR<M>::MinorDiagR(M mat, ssize_t shift) : mat(std::forward<M>(mat)), shift(shift) {
        assert(mat.isSquare());
        assert(std::abs(shift) < mat.getRow());
    }

    template<Matrix M>
    decltype(auto) MinorDiagR<M>::calc(size_t index) const noexcept {
        size_t r = shift < 0 ? -shift : 0;
        size_t c = shift > 0 ? shift : 0;
        return mat.calc(r + index, c + index);
    }

    template<Matrix M>
    auto&& MinorDiagR<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<MinorDiagR<M>> {
        using Expr = std::remove_cvref<M>::type;
    public:
        using ScalarType = Expr::ScalarType;
        constexpr static bool FastAssign = false;
    };
}
