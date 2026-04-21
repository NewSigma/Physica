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

#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Scalar T, size_t Order = Dynamic>
    class DiagMatrix : public RValueMatrix<DiagMatrix<T, Order>> {
        using This = DiagMatrix<T, Order>;
        using Base = RValueMatrix<This>;
        using VectorType = DenseVector<T, Order>;
    protected:
        using typename Base::Tr;
        using typename Base::Trv;
        using typename Base::Tv;
    private:
        VectorType diags;
    public:
        DiagMatrix() = default;
        explicit DiagMatrix(size_t order);
        explicit DiagMatrix(const Vector auto& diags_);
        DiagMatrix(const This&) = default;
        DiagMatrix(This&&) noexcept = default;
        ~DiagMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        using Base::operator*;
        void operator*=(const Scalar auto& x);
        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] Tr lnAbsDet() const;
        [[nodiscard]] T sgndet() const;

        [[nodiscard]] This inv() const;
        [[nodiscard]] auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) hermite(this auto&&) noexcept;

        void resize(size_t order);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& diag(this auto&& self) noexcept;
        [[nodiscard]] size_t getOrder() const noexcept;
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        /* Static members */
        [[nodiscard]] static This identity(size_t order);
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept { return true; }
    };

    template<Scalar T, size_t Order>
    DiagMatrix<T, Order>::DiagMatrix(size_t order) : diags(order) {}

    template<Scalar T, size_t Order>
    DiagMatrix<T, Order>::DiagMatrix(const Vector auto& diags_) : diags(diags_) {}

    template<Scalar T, size_t Order>
    void DiagMatrix<T, Order>::operator*=(const Scalar auto& x) {
        diags *= x;
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::operator*(this auto&& self, Vector auto&& v) noexcept {
        assert(self.getOrder() == v.getLength());
        return hadamard(std::forward<decltype(self)>(self).diag(), std::forward<decltype(v)>(v));
    }

    template<Scalar T, size_t Order>
    T DiagMatrix<T, Order>::calc(size_t row, size_t col) const {
        return row == col ? diags[row] : T(0);
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::calc_value(size_t row, size_t col) const -> Tv {
        return calc(row, col).value();
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::lnAbsDet() const -> Tr {
        return ln(abs(diags)).sum();
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::sgndet() const -> T {
        return unit(diags).prod();
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::inv() const -> This {
        return This(reciprocal(diags));
    }

    template<Scalar T, size_t Order>
    auto&& DiagMatrix<T, Order>::transpose(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self);
    }

    template<Scalar T, size_t Order>
    decltype(auto) DiagMatrix<T, Order>::hermite(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).conjugate();
    }

    template<Scalar T, size_t Order>
    void DiagMatrix<T, Order>::resize(size_t order) {
        diags.resize(order);
    }

    template<Scalar T, size_t Order>
    void DiagMatrix<T, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        diags.swap(obj.diags);
    }

    template<Scalar T, size_t Order>
    auto&& DiagMatrix<T, Order>::diag(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).diags;
    }

    template<Scalar T, size_t Order>
    size_t DiagMatrix<T, Order>::getOrder() const noexcept {
        if constexpr (Order == Dynamic)
            return diags.getLength();
        return Order;
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::identity(size_t order) -> This {
        return This(VectorType(order, 1));
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DiagMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Major = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#include "DiagMatrixImpl/GEMM.h"
#include "DiagMatrixImpl/InverseGEMM.h"
