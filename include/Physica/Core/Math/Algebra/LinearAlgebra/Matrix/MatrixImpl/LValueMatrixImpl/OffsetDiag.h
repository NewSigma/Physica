/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    class OffsetDiag<M> : public LValueVector<OffsetDiag<M>> {
        using This = OffsetDiag<M>;
        using Base = LValueVector<This>;
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
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept { return mat.getRow() - (offset < 0 ? -offset : offset); } // FIXME: std::abs() is not constexpr, hence manual abs
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    OffsetDiag<M>::OffsetDiag(M mat, ssize_t offset) : mat(std::forward<M>(mat)), offset(offset) {
        assert(mat.isSquare());
        assert(std::abs(offset) < mat.getRow());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    auto OffsetDiag<M>::data_ptr(this auto&& self, size_t index) noexcept {
        auto offset = self.offset;
        size_t r = offset < 0 ? -offset : 0;
        size_t c = offset > 0 ? offset : 0;
        return self.mat.data_ptr(r + index, c + index);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    auto&& OffsetDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}
