/*
 * Copyright 2022-2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType, size_t Order, size_t MaxOrder> class DenseSymmMatrix;
    template<class ScalarType, size_t Order = Dynamic, size_t MaxOrder = Order> class DenseHermiteMatrix;

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

    template<class ScalarType, size_t Order, size_t MaxOrder>
    class DenseHermiteMatrix : public RValueMatrix<DenseHermiteMatrix<ScalarType, Order, MaxOrder>>
                             , private Internal::HalfDenseMatrixStorage<ScalarType, Order, MaxOrder> {
        static_assert(ScalarType::isComplex);

        using Storage = Internal::HalfDenseMatrixStorage<ScalarType, Order, MaxOrder>;
        using RealType = typename ScalarType::RealType;
    public:
        using ColMatrix = DenseHermiteMatrix<ScalarType, Order, MaxOrder>;
        using RowMatrix = DenseHermiteMatrix<ScalarType, Order, MaxOrder>;
        using RealMatrix = DenseSymmMatrix<typename ScalarType::RealType, Order, MaxOrder>;
        constexpr static bool isComplex = true;
    public:
        template<class OtherMatrix>
        DenseHermiteMatrix(const RValueMatrix<OtherMatrix>& mat);
        using Storage::Storage;
        DenseHermiteMatrix(const DenseHermiteMatrix&) = default;
        DenseHermiteMatrix(DenseHermiteMatrix&&) noexcept = default;
        ~DenseHermiteMatrix() = default;
        /* Operators */
        DenseHermiteMatrix& operator=(DenseHermiteMatrix m) noexcept;
        DenseHermiteMatrix& operator=(RealType r);
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        /* Operations */
        using Storage::resize;
        /* Getters */
        using Storage::getRow;
        using Storage::getColumn;
        /* Helpers */
        void swap(DenseHermiteMatrix& m) noexcept;
        /* Static members */
        [[nodiscard]] static DenseHermiteMatrix unitMatrix(size_t order);
        template<class MatrixType>
        [[nodiscard]] bool isHermiteMatrix(const RValueMatrix<MatrixType>& mat, double precision);
    };
    /**
     * Save the upper triangle part of \param mat
     */
    template<class ScalarType, size_t Order, size_t MaxOrder>
    template<class OtherMatrix>
    DenseHermiteMatrix<ScalarType, Order, MaxOrder>::DenseHermiteMatrix(const RValueMatrix<OtherMatrix>& mat)
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

    template<class ScalarType, size_t Order, size_t MaxOrder>
    DenseHermiteMatrix<ScalarType, Order, MaxOrder>&
    DenseHermiteMatrix<ScalarType, Order, MaxOrder>::operator=(DenseHermiteMatrix m) noexcept {
        swap(m);
        return *this;
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    DenseHermiteMatrix<ScalarType, Order, MaxOrder>&
    DenseHermiteMatrix<ScalarType, Order, MaxOrder>::operator=(RealType r) {
        const size_t element_count = (getRow() + 1) * getRow() / 2;
        for (size_t i = 0; i < element_count; ++i)
            Storage::operator[](i) = r;
        return *this;
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    ScalarType& DenseHermiteMatrix<ScalarType, Order, MaxOrder>::operator()(size_t row, size_t col) {
        assert(row <= col); //Optimize: possible to make use of this condition
        const size_t index = Storage::accessingIndex(row, col);
        return Storage::operator[](index);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    const ScalarType& DenseHermiteMatrix<ScalarType, Order, MaxOrder>::operator()(size_t row, size_t col) const {
        const size_t index = Storage::accessingIndex(row, col);
        return Storage::operator[](index);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    ScalarType DenseHermiteMatrix<ScalarType, Order, MaxOrder>::calc(size_t row, size_t col) const {
        const size_t index = Storage::accessingIndex(row, col);
        return col >= row ? Storage::operator[](index) : Storage::operator[](index).conjugate();
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    void DenseHermiteMatrix<ScalarType, Order, MaxOrder>::swap(DenseHermiteMatrix& m) noexcept {
        Storage::swap(m);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    DenseHermiteMatrix<ScalarType, Order, MaxOrder> DenseHermiteMatrix<ScalarType, Order, MaxOrder>::unitMatrix(size_t order) {
        DenseHermiteMatrix<ScalarType, Order, MaxOrder> result(order);
        result.toUnitMatrix();
        return result;
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    template<class MatrixType>
    bool DenseHermiteMatrix<ScalarType, Order, MaxOrder>::isHermiteMatrix(const RValueMatrix<MatrixType>& mat, double precision) {
        if (mat.getRow() != mat.getColumn())
            return false;
        
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = 0; j < mat.getColumn(); ++j)
                if (!scalarNear(mat.calc(i, j), mat.calc(j, i).conjugate(), precision))
                    return false;
        return true;
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    inline void swap(
            Physica::Core::DenseHermiteMatrix<ScalarType, Order, MaxOrder>& m1,
            Physica::Core::DenseHermiteMatrix<ScalarType, Order, MaxOrder>& m2) noexcept {
        m1.swap(m2);
    }
}
