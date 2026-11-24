/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class OffsetDiag : public RValueVector<OffsetDiag<M>> {
        using This = OffsetDiag<M>;
        using Base = RValueVector<This>;
    private:
        decay_rvalue_t<M> mat;
        ssize_t offset;
    public:
        OffsetDiag(M mat, ssize_t offset);
        OffsetDiag(const This&) = default;
        OffsetDiag(This&&) noexcept = default;
        ~OffsetDiag() = default;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(size_t index) const noexcept;
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept { return mat.getRow() - (offset < 0 ? -offset : offset); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M>
    OffsetDiag<M>::OffsetDiag(M mat, ssize_t offset) : mat(std::forward<M>(mat)), offset(offset) {
        assert(mat.isSquare());
        assert(static_cast<size_t>(std::abs(offset)) < mat.getRow());
    }

    template<Matrix M>
    decltype(auto) OffsetDiag<M>::calc(size_t index) const noexcept {
        size_t r = offset < 0 ? -offset : 0;
        size_t c = offset > 0 ? offset : 0;
        return mat.calc(r + index, c + index);
    }

    template<Matrix M>
    auto&& OffsetDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<OffsetDiag<M>> {
    public:
        using ScalarType = std::remove_cvref<M>::type::ScalarType;
    };
}
