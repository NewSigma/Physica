/*
 * Copyright 2021 WeiBo He.
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

#include "RValueMatrix.h"

namespace Physica::Core {
    /**
     * \class LValueMatrix is base class of matrixes that can be assigned to \class LValueMatrix
     * and other matrixes can be assigned to this class.
     * In other words, you can take the address of elements in the matrix.
     */
    template<class Derived>
    class LValueMatrix : public RValueMatrix<Derived> {
        using Base = RValueMatrix<Derived>;
    public:
        using typename Base::ScalarType;
    public:
        /* Operators */
        [[nodiscard]] ScalarType& operator()(size_t row, size_t column) { return Base::getDerived()(row, column); }
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t column) const { return Base::getDerived()(row, column); }
        /* Operations */
        ScalarType determinate() const;
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        void rowReduce(size_t r1, size_t r2, const ScalarType& factor);
        void columnReduce(size_t c1, size_t c2, size_t elementIndex);
        void columnReduce(size_t c1, size_t c2, const ScalarType& factor);
        inline void majorReduce(size_t v1, size_t v2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, const ScalarType& factor);
        void rowMulScalar(size_t r, const ScalarType& factor);
        void columnMulScalar(size_t c, const ScalarType& factor);
        inline void majorMulScalar(size_t v, const ScalarType& factor);
        inline void majorSwap(size_t v1, size_t v2);
        /* Getters */
        [[nodiscard]] ScalarType& getElementFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] const ScalarType& getElementFromMajorMinor(size_t major, size_t minor) const;
    };
}

#include "LValueMatrixImpl.h"
