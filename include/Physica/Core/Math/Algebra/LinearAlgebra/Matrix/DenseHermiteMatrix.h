/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Utils/Container/SymmArray.h"
#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Scalar T, size_t Order> class DenseSymmMatrix;
    /**
     * \class DenseHermiteMatrix: Save the upper triangle part of hermite matrix
     */
    template<Scalar T, size_t Order = Dynamic>
    class DenseHermiteMatrix : public RValueMatrix<DenseHermiteMatrix<T, Order>> {
        using This = DenseHermiteMatrix<T, Order>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::Tr;
    private:
        SymmArray<T, Order> storage;
    public:
        DenseHermiteMatrix() = default;
        DenseHermiteMatrix(size_t order);
        DenseHermiteMatrix(size_t order, const T& t);
        DenseHermiteMatrix(const Matrix auto& mat);
        DenseHermiteMatrix(const This&) = default;
        DenseHermiteMatrix(This&&) noexcept = default;
        ~DenseHermiteMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        This& operator=(Tr value);
        [[nodiscard]] decltype(auto) operator[](size_t row, size_t col);
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] const This& hermite() const noexcept { return *this; }

        void resize(size_t order, auto&&... args);
        void swap(This& __restrict m) noexcept;

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);

        const H5DataSet<1> read(const H5Loc& loc, const char* name);
        H5DataSet<1> write(H5Loc& loc, const char* name) const;
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return storage.getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] auto&& asVector(this auto&& self) noexcept { return self.storage.asArray(); }
        /* Static members */
        [[nodiscard]] static This identity(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t order, auto& distribution);
    };

    template<Scalar T, size_t Order>
    DenseHermiteMatrix<T, Order>::DenseHermiteMatrix(size_t order) : storage(order) {}

    template<Scalar T, size_t Order>
    DenseHermiteMatrix<T, Order>::DenseHermiteMatrix(size_t order, const T& t) : storage(order, t) {}

    template<Scalar T, size_t Order>
    DenseHermiteMatrix<T, Order>::DenseHermiteMatrix(const Matrix auto& mat)
            : DenseHermiteMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getCol());
        size_t index = 0;
        for (size_t i = 0; i < mat.getRow(); ++i) {
            for (size_t j = i; j < mat.getRow(); ++j) {
                asVector()[index] = mat.calc(i, j);
                ++index;
            }
        }
    }

    template<Scalar T, size_t Order>
    auto DenseHermiteMatrix<T, Order>::operator=(Tr value) -> This& {
        asVector() = value;
        return *this;
    }

    template<Scalar T, size_t Order>
    decltype(auto) DenseHermiteMatrix<T, Order>::operator[](size_t row, size_t col) {
        assert(row <= col); // Optimize: possible to make use of this condition
        return storage[row, col];
    }

    template<Scalar T, size_t Order>
    auto DenseHermiteMatrix<T, Order>::calc(size_t row, size_t col) const -> T {
        T x = storage[row, col];
        return col >= row ? x : x.conjugate();
    }

    template<Scalar T, size_t Order>
    void DenseHermiteMatrix<T, Order>::resize(size_t order, auto&&... args) {
        storage.resize(order, std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, size_t Order>
    void DenseHermiteMatrix<T, Order>::swap(DenseHermiteMatrix& __restrict m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        storage.swap(m.storage);
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    void DenseHermiteMatrix<T, Order>::random_uniform() {
        asVector().template random_uniform<R>();
        for (size_t i = 0; i < getRow(); ++i)
            (*this)[i, i].imag() = Tr(0);
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    void DenseHermiteMatrix<T, Order>::random_normal() {
        asVector().template random_normal<R>();
        for (size_t i = 0; i < getRow(); ++i)
            (*this)[i, i].imag() = Tr(0);
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    void DenseHermiteMatrix<T, Order>::random_any(auto& distribution) {
        asVector().template random_any<R>(distribution);
        for (size_t i = 0; i < getRow(); ++i)
            (*this)[i, i].imag() = Tr(0);
    }

    template<Scalar T, size_t Order>
    auto DenseHermiteMatrix<T, Order>::identity(size_t order) -> This {
        This result(order);
        result.toIdentity();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseHermiteMatrix<T, Order>::random_uniform(size_t order) -> This {
        This result(order);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseHermiteMatrix<T, Order>::random_normal(size_t order) -> This {
        This result(order);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    auto DenseHermiteMatrix<T, Order>::random_any(size_t order, auto& distribution) -> This {
        This result(order);
        result.template random_any<R>(distribution);
        return result;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T, size_t Order>
    const H5DataSet<1> DenseHermiteMatrix<T, Order>::read(
            const H5Loc& loc, const char* name) {
        return asVector().read(loc, name);
    }

    template<Scalar T, size_t Order>
    H5DataSet<1> DenseHermiteMatrix<T, Order>::write(
            H5Loc& loc, const char* name) const {
        return asVector().write(loc, name);
    }
#endif

    template<Scalar T, size_t Order>
    void swap(DenseHermiteMatrix<T, Order>& __restrict m1, DenseHermiteMatrix<T, Order>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DenseHermiteMatrix<T, Order>> {
        static_assert(T::isComplex, "[Error]: Using a symmetric matrix is preferred for real numbers");
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
