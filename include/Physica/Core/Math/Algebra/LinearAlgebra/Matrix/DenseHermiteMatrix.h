/*
 * Copyright 2022-2025 Weibo He.
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

    template<Scalar T, size_t Order = Dynamic>
    class DenseHermiteMatrix : public RValueMatrix<DenseHermiteMatrix<T, Order>>
                             , private SymmArray<T, Order> {
        using This = DenseHermiteMatrix<T, Order>;
        using Base = RValueMatrix<This>;
        using Storage = SymmArray<T, Order>;
        using VectorStorage = Storage::ArrayType;
    public:
        using typename Base::ScalarType;
        using RealType = ScalarType::RealType;
        using ColMatrix = This;
        using RowMatrix = This;
        constexpr static bool isComplex = true;
    public:
        using Storage::Storage;
        template<Matrix M>
        DenseHermiteMatrix(const M& mat);
        DenseHermiteMatrix(const This&) = default;
        DenseHermiteMatrix(This&&) noexcept = default;
        ~DenseHermiteMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        inline This& operator=(RealType value);
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        using Base::format;
        using Base::transpose;
        using Storage::resize;
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        void swap(This& __restrict m) noexcept;

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R, class Distribution>
        void random_any(Distribution& dist);

        inline const H5DataSet<1> read(const H5Loc& loc, const char* name);
        inline H5DataSet<1> write(H5Loc& loc, const char* name) const;
        /* Getters */
        using Base::getDerived;
        using Storage::getCol;
        using Storage::getRow;
        [[nodiscard]] VectorStorage& asVector() noexcept { return Storage::asArray(); }
        [[nodiscard]] const VectorStorage& asVector() const noexcept { return Storage::asArray(); }
        /* Static members */
        [[nodiscard]] static This unitMatrix(size_t order);
        template<RNG R>
        [[nodiscard]] inline static This random_uniform(size_t order);
        template<RNG R>
        [[nodiscard]] inline static This random_normal(size_t order);
        template<RNG R, class Distribution>
        [[nodiscard]] inline static This random_any(size_t order, Distribution& dist);
        template<Matrix M>
        [[nodiscard]] bool isHermiteMatrix(const M& mat, double precision);
    };
    /**
     * Save the upper triangle part of \param mat
     */
    template<Scalar T, size_t Order>
    template<Matrix M>
    DenseHermiteMatrix<T, Order>::DenseHermiteMatrix(const M& mat)
            : DenseHermiteMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getCol());
        size_t index = 0;
        for (size_t i = 0; i < mat.getRow(); ++i) {
            for (size_t j = i; j < mat.getRow(); ++j) {
                Storage::operator[](index) = mat.calc(i, j);
                ++index;
            }
        }
    }

    template<Scalar T, size_t Order>
    inline DenseHermiteMatrix<T, Order>& DenseHermiteMatrix<T, Order>::operator=(RealType value) {
        asVector() = value;
        return *this;
    }

    template<Scalar T, size_t Order>
    DenseHermiteMatrix<T, Order>::ScalarType&
    DenseHermiteMatrix<T, Order>::operator()(size_t row, size_t col) {
        assert(row <= col); // Optimize: possible to make use of this condition
        const size_t index = Storage::toIndex1D(row, col);
        return Storage::operator[](index);
    }

    template<Scalar T, size_t Order>
    const DenseHermiteMatrix<T, Order>::ScalarType&
    DenseHermiteMatrix<T, Order>::operator()(size_t row, size_t col) const {
        const size_t index = Storage::toIndex1D(row, col);
        return Storage::operator[](index);
    }

    template<Scalar T, size_t Order>
    DenseHermiteMatrix<T, Order>::ScalarType
    DenseHermiteMatrix<T, Order>::calc(size_t row, size_t col) const {
        const size_t index = Storage::toIndex1D(row, col);
        return col >= row ? Storage::operator[](index) : Storage::operator[](index).conjugate();
    }

    template<Scalar T, size_t Order>
    void DenseHermiteMatrix<T, Order>::swap(DenseHermiteMatrix& __restrict m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        Storage::swap(m);
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    void DenseHermiteMatrix<T, Order>::random_uniform() {
        asVector().template random_uniform<R>();
        for (size_t i = 0; i < getRow(); ++i)
            this->operator()(i, i).imag() = RealType(0);
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    void DenseHermiteMatrix<T, Order>::random_normal() {
        asVector().template random_normal<R>();
        for (size_t i = 0; i < getRow(); ++i)
            this->operator()(i, i).imag() = RealType(0);
    }

    template<Scalar T, size_t Order>
    template<RNG R, class Distribution>
    void DenseHermiteMatrix<T, Order>::random_any(Distribution& dist) {
        asVector().template random_any<R, Distribution>(dist);
        for (size_t i = 0; i < getRow(); ++i)
            this->operator()(i, i).imag() = RealType(0);
    }

    template<Scalar T, size_t Order>
    DenseHermiteMatrix<T, Order> DenseHermiteMatrix<T, Order>::unitMatrix(size_t order) {
        DenseHermiteMatrix<T, Order> result(order);
        result.toUnitMatrix();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    [[nodiscard]] inline DenseHermiteMatrix<T, Order>
    DenseHermiteMatrix<T, Order>::random_uniform(size_t order) {
        This result(order);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R>
    [[nodiscard]] inline DenseHermiteMatrix<T, Order>
    DenseHermiteMatrix<T, Order>::random_normal(size_t order) {
        This result(order);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Order>
    template<RNG R, class Distribution>
    [[nodiscard]] inline DenseHermiteMatrix<T, Order>
    DenseHermiteMatrix<T, Order>::random_any(size_t order, Distribution& dist) {
        This result(order);
        result.template random_any<R, Distribution>(dist);
        return result;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T, size_t Order>
    inline const H5DataSet<1> DenseHermiteMatrix<T, Order>::read(
            const H5Loc& loc, const char* name) {
        return asVector().read(loc, name);
    }

    template<Scalar T, size_t Order>
    inline H5DataSet<1> DenseHermiteMatrix<T, Order>::write(
            H5Loc& loc, const char* name) const {
        return asVector().write(loc, name);
    }
#endif
    template<Scalar T, size_t Order>
    template<Matrix M>
    bool DenseHermiteMatrix<T, Order>::isHermiteMatrix(const M& mat, double precision) {
        if (mat.getRow() != mat.getCol())
            return false;

        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = 0; j < mat.getCol(); ++j)
                if (!scalarNear(mat.calc(i, j), mat.calc(j, i).conjugate(), precision))
                    return false;
        return true;
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DenseHermiteMatrix<T, Order>> {
        static_assert(T::isComplex, "[Error]: Using a symmetric matrix is preferred for real numbers");
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

namespace std {
    template<Physica::Scalar T, size_t Order>
    inline void swap(
            Physica::DenseHermiteMatrix<T, Order>& __restrict m1,
            Physica::DenseHermiteMatrix<T, Order>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}
