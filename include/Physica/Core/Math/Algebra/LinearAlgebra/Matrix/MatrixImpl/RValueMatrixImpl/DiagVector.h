/*
 * Copyright 2025-2026 Weibo He.
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
    class DiagVectorR : public RValueVector<DiagVectorR<M>> {
        using This = DiagVectorR<M>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M> mat;
    public:
        explicit DiagVectorR(M&& mat) : mat(std::forward<M>(mat)) {}
        DiagVectorR(const This&) = default;
        DiagVectorR(This&&) = default;
        ~DiagVectorR() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const { return mat.calc(index, index); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
    };

    template<Matrix M>
    auto&& DiagVectorR<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<DiagVectorR<M>> {
        using Expr = std::remove_cvref<M>::type;
    public:
        using ScalarType = Expr::ScalarType;
        constexpr static size_t SizeAtCompile = std::max(Expr::RowAtCompile, Expr::ColAtCompile);
        constexpr static bool FastAssign = false;
    };
}
