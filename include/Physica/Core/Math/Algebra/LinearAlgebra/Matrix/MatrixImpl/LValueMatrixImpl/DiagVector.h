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

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    class DiagVector<M> : public LValueVector<DiagVector<M>> {
        using This = DiagVector<M>;
        using Base = LValueVector<This>;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit DiagVector(M&& mat) : mat(std::forward<M>(mat)) {}
        DiagVector(const This&) = default;
        DiagVector(This&&) = default;
        ~DiagVector() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept { return self.mat.data_ptr(index, index); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    auto&& DiagVector<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    size_t DiagVector<M>::getLength() const noexcept {
        return std::min(mat.getCol(), mat.getRow());
    }

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    __host__ __device__ consteval size_t DiagVector<M>::getSizeAtCompile() noexcept {
        using Expr = std::remove_cvref<M>::type;
        return std::max(Expr::getRowAtCompile(), Expr::getColAtCompile());
    }
}
