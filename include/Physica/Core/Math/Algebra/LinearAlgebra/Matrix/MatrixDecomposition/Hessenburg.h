/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h"

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
     * Decomposite matrix A like A = QHQ^H
     * 
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class MatrixType>
    class Hessenburg {
        static_assert(MatrixType::RowAtCompile > 2 || MatrixType::RowAtCompile == Dynamic,
                      "Unnecessary hessenburg operation on matrixes whose rank is 1 or 2");
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;
        using MatrixH = HessenburgMatrixH<MatrixType>;
        using WorkingMatrix = typename MatrixType::ColMatrix;

        constexpr static size_t normVectorLength = MatrixType::RowAtCompile == Dynamic ? Dynamic : (MatrixType::RowAtCompile - 2);
        using HouseholderNorm = Vector<ScalarType, normVectorLength, normVectorLength>;
    private:
        WorkingMatrix working;
        HouseholderNorm normVector;
        const MatrixType& source;
    public:
        Hessenburg(const LValueMatrix<MatrixType>& source_);
        /* Getters */
        [[nodiscard]] MatrixH getMatrixH() const noexcept { return MatrixH(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixQ() const noexcept;
    private:
        void compute();
        friend class HessenburgMatrixH<MatrixType>;
    };

    template<class MatrixType>
    Hessenburg<MatrixType>::Hessenburg(const LValueMatrix<MatrixType>& source_)
            : working(source_.getRow(), source_.getRow())
            , normVector(source_.getRow() - 2)
            , source(source_.getDerived()) {
        assert(source.getRow() == source.getColumn());
        assert(source.getRow() > 2);
        compute();
    }

    template<class MatrixType>
    HouseholderSequence<typename Hessenburg<MatrixType>::WorkingMatrix> Hessenburg<MatrixType>::getMatrixQ() const noexcept {
        HouseholderSequence result(working);
        result.setSize(working.getRow() - 2);
        result.setShift(1);
        return result;
    }

    template<class MatrixType>
    void Hessenburg<MatrixType>::compute() {
        const size_t order = source.getRow();
        working = source;
        for (size_t i = 0; i < order - 2; ++i) {
            auto to_col = working.col(i);
            auto temp = to_col.tail(i + 1);

            ScalarType unit;
            if (temp[0].squaredNorm() <= RealType(std::numeric_limits<RealType>::min()))
                unit = ScalarType::One();
            else
                unit = temp[0].unit();

            const RealType norm = householderInPlace(temp);
            normVector[i] = -norm * unit;

            if (!norm.isZero()) {
                auto target_right = working.rightCols(i + 1);
                applyHouseholder(target_right, temp);
                auto target_bottomRight = working.bottomRightCorner(i + 1);
                applyHouseholder(temp, target_bottomRight);
            }
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
            auto copy = toCol.head(i + 1);
            copy = fromCol.head(i + 1);
            toCol[i + 1] = hess.normVector[i];
            auto zero = toCol.tail(i + 2);
            zero = ScalarType::Zero();
        }
        auto copy = target.rightCols(i);
        copy = hess.working.rightCols(i);
    }
}
