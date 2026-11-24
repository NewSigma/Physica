/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/StridedVector.h"
#include "../StridedMatrix.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    class MainDiag<M> : public StridedVector<MainDiag<M>> {
        using This = MainDiag<M>;
        using Base = StridedVector<This>;
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
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] size_t getStride() const noexcept;
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    auto&& MainDiag<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    constexpr size_t MainDiag<M>::getLength() const noexcept {
        return std::min(mat.getCol(), mat.getRow());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    auto MainDiag<M>::data_handle(this auto&& self) noexcept {
        return self.mat.data_handle();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    size_t MainDiag<M>::getStride() const noexcept {
        return mat.getRowStride() + mat.getColStride();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ consteval size_t MainDiag<M>::getSizeAtCompile() noexcept {
        using Expr = std::remove_cvref<M>::type;
        return std::min(Expr::getRowAtCompile(), Expr::getColAtCompile());
    }
}
