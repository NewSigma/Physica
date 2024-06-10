/*
 * Copyright 2022-2024 WeiBo He.
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
#include "DenseMatrixImpl/HalfDenseMatrixStorage.h"
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixProduct/MatrixProduct.h>

namespace Physica::Core {
    template<class T, size_t Order, size_t MaxOrder> class DenseSymmMatrix;
    template<class T, size_t Order = Dynamic, size_t MaxOrder = Order> class DenseHermiteMatrix;

    namespace Internal {
        template<class T> class Traits;

        template<class T, size_t Order, size_t MaxOrder>
        class Traits<DenseHermiteMatrix<T, Order, MaxOrder>> {
        public:
            using ScalarType = T;
            constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
            constexpr static size_t RowAtCompile = Order;
            constexpr static size_t ColumnAtCompile = Order;
            constexpr static size_t MaxRowAtCompile = MaxOrder;
            constexpr static size_t MaxColumnAtCompile = MaxOrder;
            constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
        };
    }

    template<class T, size_t Order, size_t MaxOrder>
    class DenseHermiteMatrix : public RValueMatrix<DenseHermiteMatrix<T, Order, MaxOrder>>
                             , private Internal::HalfDenseMatrixStorage<T, Order, MaxOrder> {
        static_assert(T::isComplex, "[Error]: Using a symmetric matrix is preferred for real numbers");
        using This = DenseHermiteMatrix<T, Order, MaxOrder>;
        using Base = RValueMatrix<This>;
        using Storage = Internal::HalfDenseMatrixStorage<T, Order, MaxOrder>;
        using VectorBase = typename Storage::Base;
    public:
        using typename Base::ScalarType;
        using RealType = typename ScalarType::RealType;
        using ColMatrix = This;
        using RowMatrix = This;
        using RealMatrix = DenseSymmMatrix<typename ScalarType::RealType, Order, MaxOrder>;
        constexpr static bool isComplex = true;
    public:
        template<class OtherMatrix>
        DenseHermiteMatrix(const RValueMatrix<OtherMatrix>& mat);
        using Storage::Storage;
        DenseHermiteMatrix(const This&) = default;
        DenseHermiteMatrix(This&&) noexcept = default;
        ~DenseHermiteMatrix() = default;
        /* Operators */
        DenseHermiteMatrix& operator=(This m) noexcept;
        DenseHermiteMatrix& operator=(RealType r);
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
        template<class VectorType>
        [[nodiscard]] inline MatrixVectorProduct<This, VectorType> operator*(const RValueVector<VectorType>& vec) const noexcept;
        /* Operations */
        using Base::transpose;
        using Storage::resize;
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        void swap(This& __restrict m) noexcept;

        template<class RandomGenerator>
        void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);
        /* Getters */
        using Base::getDerived;
        using Storage::getRow;
        using Storage::getColumn;
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
        /* Static members */
        [[nodiscard]] static This unitMatrix(size_t order);
        template<class MatrixType>
        [[nodiscard]] bool isHermiteMatrix(const RValueMatrix<MatrixType>& mat, double precision);
    };
    /**
     * Save the upper triangle part of \param mat
     */
    template<class T, size_t Order, size_t MaxOrder>
    template<class OtherMatrix>
    DenseHermiteMatrix<T, Order, MaxOrder>::DenseHermiteMatrix(const RValueMatrix<OtherMatrix>& mat)
            : DenseHermiteMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getColumn());
        size_t index = 0;
        for (size_t i = 0; i < mat.getRow(); ++i) {
            for (size_t j = i; j < mat.getRow(); ++j) {
                Storage::operator[](index) = mat.calc(i, j);
                ++index;
            }
        }
    }

    template<class T, size_t Order, size_t MaxOrder>
    DenseHermiteMatrix<T, Order, MaxOrder>&
    DenseHermiteMatrix<T, Order, MaxOrder>::operator=(DenseHermiteMatrix m) noexcept {
        swap(m);
        return *this;
    }

    template<class T, size_t Order, size_t MaxOrder>
    DenseHermiteMatrix<T, Order, MaxOrder>&
    DenseHermiteMatrix<T, Order, MaxOrder>::operator=(RealType r) {
        const size_t element_count = (getRow() + 1) * getRow() / 2;
        for (size_t i = 0; i < element_count; ++i)
            Storage::operator[](i) = r;
        return *this;
    }

    template<class T, size_t Order, size_t MaxOrder>
    typename DenseHermiteMatrix<T, Order, MaxOrder>::ScalarType&
    DenseHermiteMatrix<T, Order, MaxOrder>::operator()(size_t row, size_t col) {
        assert(row <= col); //Optimize: possible to make use of this condition
        const size_t index = Storage::toIndex1D(row, col);
        return Storage::operator[](index);
    }

    template<class T, size_t Order, size_t MaxOrder>
    const typename DenseHermiteMatrix<T, Order, MaxOrder>::ScalarType&
    DenseHermiteMatrix<T, Order, MaxOrder>::operator()(size_t row, size_t col) const {
        const size_t index = Storage::toIndex1D(row, col);
        return Storage::operator[](index);
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class VectorType>
    inline MatrixVectorProduct<DenseHermiteMatrix<T, Order, MaxOrder>, VectorType>
    DenseHermiteMatrix<T, Order, MaxOrder>::operator*(const RValueVector<VectorType>& vec) const noexcept {
        return static_cast<const Base&>(*this) * vec;
    }

    template<class T, size_t Order, size_t MaxOrder>
    typename DenseHermiteMatrix<T, Order, MaxOrder>::ScalarType
    DenseHermiteMatrix<T, Order, MaxOrder>::calc(size_t row, size_t col) const {
        const size_t index = Storage::toIndex1D(row, col);
        return col >= row ? Storage::operator[](index) : Storage::operator[](index).conjugate();
    }

    template<class T, size_t Order, size_t MaxOrder>
    void DenseHermiteMatrix<T, Order, MaxOrder>::swap(DenseHermiteMatrix& __restrict m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        Storage::swap(m);
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class RandomGenerator>
    void DenseHermiteMatrix<T, Order, MaxOrder>::random_uniform(RandomGenerator& gen) {
        Storage::random_uniform(gen);
        for (size_t i = 0; i < getRow(); ++i)
            this->operator()(i, i).getImag() = RealType(0);
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class RandomGenerator>
    void DenseHermiteMatrix<T, Order, MaxOrder>::random_normal(RandomGenerator& gen) {
        Storage::random_normal(gen);
        for (size_t i = 0; i < getRow(); ++i)
            this->operator()(i, i).getImag() = RealType(0);
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class Distribution, class RandomGenerator>
    void DenseHermiteMatrix<T, Order, MaxOrder>::random_any(Distribution& dist, RandomGenerator& gen) {
        Storage::random_any(dist, gen);
        for (size_t i = 0; i < getRow(); ++i)
            this->operator()(i, i).getImag() = RealType(0);
    }

    template<class T, size_t Order, size_t MaxOrder>
    DenseHermiteMatrix<T, Order, MaxOrder> DenseHermiteMatrix<T, Order, MaxOrder>::unitMatrix(size_t order) {
        DenseHermiteMatrix<T, Order, MaxOrder> result(order);
        result.toUnitMatrix();
        return result;
    }

    template<class T, size_t Order, size_t MaxOrder>
    template<class MatrixType>
    bool DenseHermiteMatrix<T, Order, MaxOrder>::isHermiteMatrix(const RValueMatrix<MatrixType>& mat, double precision) {
        if (mat.getRow() != mat.getColumn())
            return false;
        
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = 0; j < mat.getColumn(); ++j)
                if (!scalarNear(mat.calc(i, j), mat.calc(j, i).conjugate(), precision))
                    return false;
        return true;
    }
}

namespace std {
    template<class T, size_t Order, size_t MaxOrder>
    inline void swap(
            Physica::Core::DenseHermiteMatrix<T, Order, MaxOrder>& __restrict m1,
            Physica::Core::DenseHermiteMatrix<T, Order, MaxOrder>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}
