/*
 * Copyright 2020-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"

namespace Physica::Core {
    /**
     * Construct using a lower triangular matrix, echo column represents a householder transformation
     *
     * \tparam ColWiseRead
     * Read data from columns of source, otherwise read from rows
     *
     * Reference:
     * [1] Eigen https://eigen.tuxfamily.org/
     */
    template<Matrix T, bool ColWiseRead = true>
    class HouseholderSequence : public RValueMatrix<HouseholderSequence<T, ColWiseRead>> {
        const T& source;
        /**
         * Number of householder transformation in this sequence
         */
        size_t size;
        /**
         * The start index of each householder vector
         */
        size_t shift;
    public:
        HouseholderSequence(const T& source_);
        /* Operations */
        template<Matrix U>
        void assignTo(LValueMatrix<U>& target) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return ColWiseRead ? source.getRow() : source.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getRow(); }
        [[nodiscard]] size_t getSize() const noexcept { return size; }
        [[nodiscard]] size_t getShift() const noexcept { return shift; }
        /* Setters */
        inline void setSize(size_t size_);
        inline void setShift(size_t shift_);
    };

    template<Matrix T, bool ColWiseRead>
    HouseholderSequence<T, ColWiseRead>::HouseholderSequence(const T& source_)
            : source(source_)
            , size(source_.getCol())
            , shift(0) {}

    template<Matrix T, bool ColWiseRead>
    template<Matrix U>
    void HouseholderSequence<T, ColWiseRead>::assignTo(LValueMatrix<U>& target) const {
        target.toUnitMatrix();
        const size_t shift1 = shift + target.getRow() - (ColWiseRead ? source.getRow() : source.getCol());
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

    template<Matrix T, bool ColWiseRead>
    inline void HouseholderSequence<T, ColWiseRead>::setSize(size_t size_) {
        assert(size_ <= source.getCol());
        size = size_;
    }

    template<Matrix T, bool ColWiseRead>
    inline void HouseholderSequence<T, ColWiseRead>::setShift(size_t shift_) {
        shift = shift_;
    }
}

namespace Physica {
    template<Matrix T, bool ColWiseRead>
    class Traits<Core::HouseholderSequence<T, ColWiseRead>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;
    };
}
