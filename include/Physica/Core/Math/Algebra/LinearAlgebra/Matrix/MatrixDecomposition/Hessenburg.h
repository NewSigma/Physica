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
    template<class MatrixType> class HessenburgMatrixH;

    namespace Internal {
        template<class MatrixType>
        class Traits<HessenburgMatrixH<MatrixType>> : public Traits<MatrixType> {
        private:
            using Base = Traits<MatrixType>;
            using Base::MatrixOption;
        };
    }
    /**
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class MatrixType>
    class Hessenburg {
        using ScalarType = typename MatrixType::ScalarType;
        using MatrixH = HessenburgMatrixH<MatrixType>;
        using WorkingMatrix = typename MatrixType::ColMatrix;
        const MatrixType& source;
        WorkingMatrix working;
    public:
        Hessenburg(const LValueMatrix<MatrixType>& source_);
        /* Getters */
        [[nodiscard]] MatrixH getMatrixH() const noexcept { return MatrixH(*this); }
    private:
        void compute();
        friend class HessenburgMatrixH<MatrixType>;
    };

    template<class MatrixType>
    Hessenburg<MatrixType>::Hessenburg(const LValueMatrix<MatrixType>& source_)
            : source(source_.getDerived())
            , working(source_.getRow(), source_.getRow()) {
        assert(source.getRow() == source.getColumn());
        compute();
    }

    template<class MatrixType>
    void Hessenburg<MatrixType>::compute() {
        constexpr static size_t householderLength = MatrixType::RowAtCompile == Dynamic ? Dynamic : (MatrixType::RowAtCompile - 1);
        using VectorType = Vector<ScalarType, householderLength, householderLength>;
        const size_t order = source.getRow();
        working = source;
        VectorType householderVector = VectorType(order - 1);
        for (size_t i = 0; i < order - 2; ++i) {
            auto to_col = working.col(i);
            auto eliminate = to_col.tail(i + 1);
            auto temp = householderVector.head(order - i - 1);
            const auto norm = householder(eliminate, temp);
            eliminate[0] = eliminate[0].isNegative() ? norm : -norm;
            eliminate.tail(1) = ScalarType::Zero();

            auto target_right = working.rightCols(i + 1);
            applyHouseholder(target_right, temp);
            auto target_bottomRight = working.bottomRightCorner(i + 1);
            applyHouseholder(temp, target_bottomRight);
        }
    }

    template<class MatrixType>
    class HessenburgMatrixH : public RValueMatrix<HessenburgMatrixH<MatrixType>> {
        using Base = RValueMatrix<HessenburgMatrixH<MatrixType>>;
        using typename Base::ScalarType;
        const Hessenburg<MatrixType>& hess;
    public:
        HessenburgMatrixH(const Hessenburg<MatrixType>& hess_) : hess(hess_) {}
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return hess.source.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return hess.source.getColumn(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void HessenburgMatrixH<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        const size_t order = getRow();
        size_t i = 0;
        for (; i < order - 2; ++i) {
            auto fromCol = hess.working.col(i);
            auto toCol = target.col(i);
            auto copy = toCol.head(i + 2);
            copy = fromCol.head(i + 2);
            auto zero = toCol.tail(i + 2);
            zero = ScalarType::Zero();
        }
        auto copy = target.rightCols(i);
        copy = hess.working.rightCols(i);
    }
}
