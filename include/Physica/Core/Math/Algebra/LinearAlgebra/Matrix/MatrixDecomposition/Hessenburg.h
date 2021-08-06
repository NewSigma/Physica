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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/RValueMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"

namespace Physica::Core {
    template<class MatrixType> class Hessenburg;

    namespace Internal {
        template<class MatrixType>
        class Traits<Hessenburg<MatrixType>> : public Traits<MatrixType> {
        private:
            using Base = Traits<MatrixType>;
            using Base::MatrixOption;
        };
    }

    template<class MatrixType>
    class Hessenburg : public RValueMatrix<Hessenburg<MatrixType>> {
        using Base = RValueMatrix<Hessenburg<MatrixType>>;
        const MatrixType& source;
    public:
        Hessenburg(const LValueMatrix<MatrixType>& source_) : source(source_.getDerived()) { assert(source.getRow() == source.getColumn()); }
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return source.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return source.getColumn(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void Hessenburg<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        constexpr static size_t householderLength = Base::RowAtCompile == Dynamic ? Dynamic : (Base::RowAtCompile - 1);
        using VectorType = Vector<typename Base::ScalarType, householderLength, householderLength>;
        const size_t order = getRow();
        target.rightCols(1) = source.rightCols(1);
        VectorType householderVector = VectorType(order - 1);
        for (size_t i = 0; i < order - 2; ++i) {
            auto from_col = source.col(i);
            auto to_col = target.col(i);
            auto temp = householderVector.head(order - i - 1);
            to_col.head(i + 1) = from_col.head(i + 1);
            to_col[0] = householder(from_col.tail(i + 1), temp);
            to_col.tail(1) = Base::ScalarType::Zero();

            auto target_right = target.rightCols(i + 1);
            applyHouseholder(target_right, temp);
            auto target_bottomRight = target.bottomRightCorner(i + 1);
            applyHouseholder(temp, target_bottomRight);
        }
    }
}
