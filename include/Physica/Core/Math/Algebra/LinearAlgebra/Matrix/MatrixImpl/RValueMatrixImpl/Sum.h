/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../RValueMatrix.h"

namespace Physica {
    template<class M, bool ReduceCol>
    class MatrixSum : public RValueVector<MatrixSum<M, ReduceCol>> {
        using This = MatrixSum<M, ReduceCol>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    private:
        const M& mat;
    public:
        MatrixSum(const M& mat_) : mat(mat_) {}
        MatrixSum(const This&) = default;
        MatrixSum(This&&) noexcept = default;
        ~MatrixSum() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<class M, bool ReduceCol>
    auto MatrixSum<M, ReduceCol>::calc(size_t index) const -> T {
        if constexpr (ReduceCol)
            return mat.row(index).sum();
        else
            return mat.col(index).sum();
    }

    template<class M, bool ReduceCol>
    auto MatrixSum<M, ReduceCol>::calc_value(size_t index) const -> Tv {
        if constexpr (ReduceCol)
            return mat.values().row(index).sum();
        else
            return mat.values().col(index).sum();
    }

    template<class M, bool ReduceCol>
    void MatrixSum<M, ReduceCol>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        if constexpr (ReduceCol)
            mat.reverse(grad);
        else
            mat.transpose().reverse(grad);
    }

    template<class M, bool ReduceCol>
    auto MatrixSum<M, ReduceCol>::values(this auto&& self) noexcept {
        return MatrixSum<decltype(self.mat.values()), ReduceCol>(self.mat.values());
    }

    template<class M, bool ReduceCol>
    constexpr size_t MatrixSum<M, ReduceCol>::getLength() const noexcept {
        if constexpr (ReduceCol)
            return mat.getRow();
        else
            return mat.getCol();
    }

    template<class M, bool ReduceCol>
    __host__ __device__ consteval size_t MatrixSum<M, ReduceCol>::getSizeAtCompile() noexcept {
        return ReduceCol ? M::getRowAtCompile() : M::getColAtCompile();
    }
}

namespace Physica {
    template<Matrix M, bool ReduceCol>
    class Traits<MatrixSum<M, ReduceCol>> {
    public:
        using ScalarType = M::ScalarType;
    };
}
