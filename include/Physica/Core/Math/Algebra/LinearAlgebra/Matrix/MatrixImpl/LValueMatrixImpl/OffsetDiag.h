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
    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    class OffsetDiag<M> : public LValueVector<OffsetDiag<M>> {
        using This = OffsetDiag<M>;
        using Base = LValueVector<This>;
    private:
        decay_rvalue_t<M> mat;
        ssize_t offset;
    public:
        OffsetDiag(M mat, ssize_t offset);
        OffsetDiag(const This&) = default;
        OffsetDiag(This&&) = default;
        ~OffsetDiag() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() - std::abs(offset); }
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    OffsetDiag<M>::OffsetDiag(M mat, ssize_t offset) : mat(std::forward<M>(mat)), offset(offset) {
        assert(mat.isSquare());
        assert(std::abs(offset) < mat.getRow());
    }

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    auto OffsetDiag<M>::data_ptr(this auto&& self, size_t index) noexcept {
        auto offset = self.offset;
        size_t r = offset < 0 ? -offset : 0;
        size_t c = offset > 0 ? offset : 0;
        return self.mat.data_ptr(r + index, c + index);
    }

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix())
    auto&& OffsetDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}
