/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    class MainDiag<M> : public LValueVector<MainDiag<M>> {
        using This = MainDiag<M>;
        using Base = LValueVector<This>;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit MainDiag(M&& mat) : mat(std::forward<M>(mat)) {}
        MainDiag(const This&) = default;
        MainDiag(This&&) noexcept = default;
        ~MainDiag() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept { return self.mat.data_ptr(index, index); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    auto&& MainDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    constexpr size_t MainDiag<M>::getLength() const noexcept {
        return std::min(mat.getCol(), mat.getRow());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isStrided())
    __host__ __device__ consteval size_t MainDiag<M>::getSizeAtCompile() noexcept {
        using Expr = std::remove_cvref<M>::type;
        return std::min(Expr::getRowAtCompile(), Expr::getColAtCompile());
    }
}
