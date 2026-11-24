/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class MainDiag : public RValueVector<MainDiag<M>> {
        using This = MainDiag<M>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit MainDiag(M&& mat) : mat(std::forward<M>(mat)) {}
        MainDiag(const This&) = default;
        MainDiag(This&&) noexcept = default;
        ~MainDiag() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const { return mat.calc(index, index); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M>
    auto&& MainDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M>
    constexpr size_t MainDiag<M>::getLength() const noexcept {
        return std::min(mat.getCol(), mat.getRow());
    }

    template<Matrix M>
    __host__ __device__ consteval size_t MainDiag<M>::getSizeAtCompile() noexcept {
        using Expr = std::remove_cvref<M>::type;
        return std::min(Expr::getRowAtCompile(), Expr::getColAtCompile());
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<MainDiag<M>> {
    public:
        using ScalarType = std::remove_cvref<M>::type::ScalarType;
    };
}
