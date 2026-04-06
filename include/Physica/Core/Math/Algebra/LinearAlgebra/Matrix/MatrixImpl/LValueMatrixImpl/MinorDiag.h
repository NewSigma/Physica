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
    class MinorDiagL final : public LValueVector<MinorDiagL<M>> {
        using This = MinorDiagL<M>;
        using Base = LValueVector<This>;
    private:
        LazyDestroy<M> mat;
        ssize_t shift;
    public:
        MinorDiagL(M mat, ssize_t shift);
        MinorDiagL(const This&) = default;
        MinorDiagL(This&&) = default;
        ~MinorDiagL() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() - std::abs(shift); }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
    };

    template<Matrix M>
    MinorDiagL<M>::MinorDiagL(M mat, ssize_t shift) : mat(std::forward<M>(mat)), shift(shift) {
        assert(mat.isSquare());
        assert(std::abs(shift) < mat.getRow());
    }

    template<Matrix M>
    auto MinorDiagL<M>::data_ptr(this auto&& self, size_t index) noexcept {
        auto shift = self.shift;
        size_t r = shift < 0 ? -shift : 0;
        size_t c = shift > 0 ? shift : 0;
        return self.mat.data_ptr(r + index, c + index);
    }

    template<Matrix M>
    auto&& MinorDiagL<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<MinorDiagL<M>> {
        using Expr = std::remove_cvref<M>::type;
    public:
        using ScalarType = Expr::ScalarType;
        constexpr static size_t SizeAtCompile = std::max(Expr::RowAtCompile, Expr::ColAtCompile);
        constexpr static bool FastAssign = false;
    };
}
