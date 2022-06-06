/*
 * Copyright 2020-2021 WeiBo He.
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

#include "LValueMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"

namespace Physica::Core {
    template<class MatrixType, bool ColWiseRead> class HouseholderSequence;

    namespace Internal {
        template<class MatrixType, bool ColWiseRead>
        class Traits<HouseholderSequence<MatrixType, ColWiseRead>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t RowAtCompile = Dynamic;
            constexpr static size_t ColumnAtCompile = Dynamic;
            constexpr static size_t MaxRowAtCompile = MatrixType::MaxRowAtCompile;
            constexpr static size_t MaxColumnAtCompile = MaxRowAtCompile;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
        };
    }
    /**
     * Construct using a lower triangular matrix, echo column represents a householder transformation
     * 
     * \tparam ColWiseRead
     * Read data from columns of source, otherwise read from rows
     * 
     * Reference:
     * [1] Eigen https://eigen.tuxfamily.org/
     */
    template<class MatrixType, bool ColWiseRead = true>
    class HouseholderSequence : public RValueMatrix<HouseholderSequence<MatrixType, ColWiseRead>> {
        const MatrixType& source;
        /**
         * Number of householder transformation in this sequence
         */
        size_t size;
        /**
         * The start index of each householder vector
         */
        size_t shift;
    public:
        HouseholderSequence(const RValueMatrix<MatrixType>& source_);
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return ColWiseRead ? source.getRow() : source.getColumn(); }
        [[nodiscard]] size_t getColumn() const noexcept { return getRow(); }
        [[nodiscard]] size_t getSize() const noexcept { return size; }
        [[nodiscard]] size_t getShift() const noexcept { return shift; }
        /* Setters */
        inline void setSize(size_t size_);
        inline void setShift(size_t shift_);
    };

    template<class MatrixType, bool ColWiseRead>
    HouseholderSequence<MatrixType, ColWiseRead>::HouseholderSequence(const RValueMatrix<MatrixType>& source_)
            : source(source_.getDerived())
            , size(source.getColumn())
            , shift(0) {}

    template<class MatrixType, bool ColWiseRead>
    template<class OtherMatrix>
    void HouseholderSequence<MatrixType, ColWiseRead>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        target.toUnitMatrix();
        const size_t shift1 = shift + target.getRow() - (ColWiseRead ? source.getRow() : source.getColumn());
        assert(shift1 < target.getRow());

        for (size_t i = 0; i < size; ++i) {
            auto block = target.rightCols(i + shift1);
            if constexpr (ColWiseRead) {
                auto col = source.col(i);
                applyHouseholder(block, col.tail(i + shift));
            }
            else {
                auto row = source.row(i);
                applyHouseholder(block, row.tail(i + shift));
            }
        }
    }

    template<class MatrixType, bool ColWiseRead>
    inline void HouseholderSequence<MatrixType, ColWiseRead>::setSize(size_t size_) {
        assert(size_ <= source.getColumn());
        size = size_;
    }

    template<class MatrixType, bool ColWiseRead>
    inline void HouseholderSequence<MatrixType, ColWiseRead>::setShift(size_t shift_) {
        shift = shift_;
    }
}
