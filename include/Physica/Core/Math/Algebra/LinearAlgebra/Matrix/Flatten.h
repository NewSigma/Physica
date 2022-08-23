/*
 * Copyright 2022 WeiBo He.
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
    template<class MatrixType> class Flatten;

    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<Flatten<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile * MatrixType::ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = MatrixType::MaxRowAtCompile * MatrixType::MaxColumnAtCompile;
        };
    }

    template<class MatrixType>
    class Flatten : public RValueVector<Flatten<MatrixType>> {
        const MatrixType& mat;
    public:
        using Base = RValueVector<Flatten<MatrixType>>;
        using typename Base::ScalarType;
    public:
        Flatten(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        template<class OtherVector>
        void assignTo(LValueVector<OtherVector>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getColumn(); }
    };

    template<class VectorType>
    template<class OtherVector>
    void Flatten<VectorType>::assignTo(LValueVector<OtherVector>& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);
    }

    template<class MatrixType>
    typename Flatten<MatrixType>::ScalarType Flatten<MatrixType>::calc(size_t index) const {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor);
    }
}
