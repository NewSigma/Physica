/*
 * Copyright 2020-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"
#include "Physica/Core/Utils/Container/SymmArray.h"

namespace Physica {
    template<Scalar T, size_t Order, Vector V>
    class SyMV;

    template<Scalar T, size_t Order = Dynamic>
    class DenseSymmMatrix : public LValueMatrix<DenseSymmMatrix<T, Order>> {
        using This = DenseSymmMatrix<T, Order>;
        using Base = LValueMatrix<This>;
    private:
        SymmArray<T, Order> storage;
    public:
        DenseSymmMatrix() = default;
        DenseSymmMatrix(size_t order);
        DenseSymmMatrix(size_t order, const T& t);
        DenseSymmMatrix(const Matrix auto& mat);
        DenseSymmMatrix(const This&) = default;
        DenseSymmMatrix(This&&) noexcept = default;
        ~DenseSymmMatrix() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        auto operator=(Scalar auto x) noexcept -> This&;
        void operator+=(Scalar auto x) noexcept;
        void operator-=(Scalar auto x) noexcept;
        void operator*=(Scalar auto x) noexcept;
        void operator/=(Scalar auto x) noexcept;

        template<Vector V>
        [[nodiscard]] auto operator*(const V& v) const noexcept;
        /* Operations */
        void resize(size_t order, auto&&... args);
        void resize(const Matrix auto& m, auto&&... args);

        [[nodiscard]] T max() const { return asVector().max(); }
        [[nodiscard]] T min() const { return asVector().min(); }
        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        [[nodiscard]] decltype(auto) hermite() const noexcept;
        void swap(This& __restrict m) noexcept;

        void zeros() { storage.zeros(); }
        template<RNG R>
        void random_uniform() { asVector().template random_uniform<R>(); }
        template<RNG R>
        void random_normal() { asVector().template random_normal<R>(); }
        template<RNG R>
        void random_any(auto& distribution) { asVector().template random_any<R>(distribution); }

        const H5DataSet<1> read(const H5Loc& loc, const char* name) { return storage.read(loc, name); }
        H5DataSet<1> write(H5Loc& loc, const char* name) const { return storage.write(loc, name); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return storage.getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t toIndex1D(size_t r, size_t c) const noexcept { return storage.toIndex1D(r, c); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] auto&& asVector(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] static DenseSymmMatrix identity(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t order, auto& distribution);
    };

    template<Scalar T, size_t Order>
    DenseSymmMatrix<T, Order>::DenseSymmMatrix(size_t order) : storage(order) {}

    template<Scalar T, size_t Order>
    DenseSymmMatrix<T, Order>::DenseSymmMatrix(size_t order, const T& t) : storage(order, t) {}
    /**
     * Assuming mat is a symmetric matrix, if it is not the case, only half of the elements is saved correctly
     */
    template<Scalar T, size_t Order>
    DenseSymmMatrix<T, Order>::DenseSymmMatrix(const Matrix auto& mat) : DenseSymmMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getCol());
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = i; j < mat.getRow(); ++j)
                storage[i, j] = mat.calc(i, j);
    }

    template<Scalar T, size_t Order>
    auto DenseSymmMatrix<T, Order>::operator=(Scalar auto x) noexcept -> This& {
        asVector() = x;
        return *this;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator+=(Scalar auto x) noexcept {
        asVector() += x;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator-=(Scalar auto x) noexcept {
        asVector() -= x;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator*=(Scalar auto x) noexcept {
        asVector() *= x;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::operator/=(Scalar auto x) noexcept {
        asVector() /= x;
    }

    template<Scalar T, size_t Order>
    template<Vector V>
    auto DenseSymmMatrix<T, Order>::operator*(const V& v) const noexcept {
        return SyMV<T, Order, V>(*this, v);
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::resize(size_t order, auto&&... args) {
        storage.resize(order, std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::resize(const Matrix auto& m, auto&&... args) {
        assert(m.isSquare());
        resize(m.getRow(), std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, size_t Order>
    decltype(auto) DenseSymmMatrix<T, Order>::hermite() const noexcept {
        if constexpr (T::isComplex)
            return Base::hermite();
        else
            return *this;
    }

    template<Scalar T, size_t Order>
    void DenseSymmMatrix<T, Order>::swap(This& __restrict m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        storage.swap(m.storage);
    }

    template<Scalar T, size_t Order>
    __host__ __device__ auto DenseSymmMatrix<T, Order>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        return self.storage.data_ptr(row, col);
    }

    template<Scalar T, size_t Order>
    auto&& DenseSymmMatrix<T, Order>::asVector(this auto&& self) noexcept {
        return self.storage.asArray();
    }

    template<Scalar T, size_t Order>
    auto DenseSymmMatrix<T, Order>::identity(size_t order) -> This {
        DenseSymmMatrix<T, Order> result(order);
        result.toIdentity();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseSymmMatrix<T, Order>::random_uniform(size_t order) -> This {
        This result(order);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseSymmMatrix<T, Order>::random_normal(size_t order) -> This {
        This result(order);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseSymmMatrix<T, Order>::random_any(size_t order, auto& distribution) -> This {
        This result(order);
        result.template random_any<R>(distribution);
        return result;
    }

    template<Physica::Scalar T, size_t Order>
    void swap(DenseSymmMatrix<T, Order>& __restrict m1, DenseSymmMatrix<T, Order>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DenseSymmMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#include "MatrixImpl/MatrixProduct/SyMV.h"
