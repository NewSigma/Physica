/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/StridedVector.h"
#include "../StridedMatrix.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    class OffsetDiag<M> : public StridedVector<OffsetDiag<M>> {
        using This = OffsetDiag<M>;
        using Base = StridedVector<This>;
    private:
        decay_rvalue_t<M> mat;
        ssize_t offset;
    public:
        OffsetDiag(M mat, ssize_t offset);
        OffsetDiag(const This&) = default;
        OffsetDiag(This&&) noexcept = default;
        ~OffsetDiag() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] size_t getStride() const noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept { return mat.getRow() - (offset < 0 ? -offset : offset); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    OffsetDiag<M>::OffsetDiag(M mat, ssize_t offset) : mat(std::forward<M>(mat)), offset(offset) {
        assert(mat.isSquare());
        assert(static_cast<size_t>(std::abs(offset)) < mat.getRow());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    auto&& OffsetDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    auto OffsetDiag<M>::data_handle(this auto&& self) noexcept {
        size_t r = self.offset < 0 ? -self.offset : 0;
        size_t c = self.offset > 0 ? self.offset : 0;
        return self.mat.data_handle() + r * self.mat.getRowStride() + c * self.mat.getColStride();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    size_t OffsetDiag<M>::getStride() const noexcept {
        return mat.getRowStride() + mat.getColStride();
    }
}
