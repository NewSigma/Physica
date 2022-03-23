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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h"

namespace Physica::Core {
    template<class MatrixType> class TridiagonalMatrixT;

    namespace Internal {
        template<class MatrixType>
        class Traits<TridiagonalMatrixT<MatrixType>> : public Traits<MatrixType> {
        private:
            using Base = Traits<MatrixType>;
            using Base::MatrixOption;
        };
    }
    /**
     * Decomposite symmetric matrix A like A = QTQ^H
     * 
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.426-428
     */
    template<class MatrixType>
    class Tridiagonalization {
        static_assert(MatrixType::RowAtCompile > 2 || MatrixType::RowAtCompile == Dynamic,
                      "Unnecessary hessenburg operation on matrixes whose rank is 1 or 2");
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;
        using MatrixT = TridiagonalMatrixT<MatrixType>;
        using SymmMatrix = DenseSymmMatrix<ScalarType, MatrixType::RowAtCompile, MatrixType::MaxRowAtCompile>;
        using HermiteMatrix = DenseHermiteMatrix<ScalarType, MatrixType::RowAtCompile, MatrixType::MaxRowAtCompile>;
        using WorkingMatrix = typename std::conditional<ScalarType::isComplex, HermiteMatrix, SymmMatrix>::type;

        constexpr static size_t normVectorLength = MatrixType::RowAtCompile == Dynamic ? Dynamic : (MatrixType::RowAtCompile - 2);
        constexpr static size_t bufferLength = MatrixType::RowAtCompile == Dynamic ? Dynamic : (MatrixType::RowAtCompile - 1);
        using HouseholderNorm = Vector<ScalarType, normVectorLength, normVectorLength>;
        using BufferVector = Vector<ScalarType, bufferLength, bufferLength>;
    private:
        WorkingMatrix working;
        HouseholderNorm normBuffer;
        BufferVector buffer;
    public:
        Tridiagonalization(size_t size);
        Tridiagonalization(const LValueMatrix<MatrixType>& source);
        /* Operations */
        void compute(const LValueMatrix<MatrixType>& source);
        /* Getters */
        [[nodiscard]] MatrixT getMatrixT() const noexcept { return MatrixT(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixQ() const noexcept;

        friend class TridiagonalMatrixT<MatrixType>;
    };

    template<class MatrixType>
    Tridiagonalization<MatrixType>::Tridiagonalization(size_t size)
            : working(size, size)
            , normBuffer(size - 2)
            , buffer(size - 1) {}

    template<class MatrixType>
    Tridiagonalization<MatrixType>::Tridiagonalization(const LValueMatrix<MatrixType>& source)
            : Tridiagonalization(source.getRow()) {
        compute(source);
    }

    template<class MatrixType>
    void Tridiagonalization<MatrixType>::compute(const LValueMatrix<MatrixType>& source) {
        const size_t order = source.getRow();
        working = source;
        for (size_t i = 0; i < order - 2; ++i) {
            auto to_col = working.col(i);
            auto temp = to_col.tail(i + 1);

            ScalarType unit;
            if (temp.calc(0).squaredNorm() <= RealType(std::numeric_limits<RealType>::min()))
                unit = ScalarType::One();
            else
                unit = temp.calc(0).unit();

            auto p = buffer.head(order - i - 1);
            const RealType norm = householder(temp, p);
            for (size_t j = 0; j < p.getLength(); ++j)
                working(i, i + j + 1) = p[j].conjugate();
            normBuffer[i] = -norm * unit;

            if (!norm.isZero()) {
                const ScalarType factor = working(i, i + 1);
                working(i, i + 1) = ScalarType::One();
                auto corner = working.bottomRightCorner(i + 1);
                p = factor * (corner * temp);
                p -= (p.conjugate() * temp * factor * ScalarType(0.5)) * temp;
                for (size_t r = 0; r < corner.getRow(); ++r)
                    for (size_t c = r; c < corner.getColumn(); ++c)
                        working(r + i + 1, c + i + 1) -= temp.calc(r) * p[c].conjugate() + temp.calc(c).conjugate() * p[r];
                working(i, i + 1) = factor;
            }
        }
    }

    template<class MatrixType>
    HouseholderSequence<typename Tridiagonalization<MatrixType>::WorkingMatrix>
    Tridiagonalization<MatrixType>::getMatrixQ() const noexcept {
        HouseholderSequence result(working);
        result.setSize(working.getRow() - 1);
        result.setShift(1);
        return result;
    }

    template<class MatrixType>
    class TridiagonalMatrixT : public RValueMatrix<TridiagonalMatrixT<MatrixType>> {
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;

        const Tridiagonalization<MatrixType>& tri;
    public:
        TridiagonalMatrixT(const Tridiagonalization<MatrixType>& tri_);
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return tri.working.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return tri.working.getColumn(); }
    };

    template<class MatrixType>
    TridiagonalMatrixT<MatrixType>::TridiagonalMatrixT(const Tridiagonalization<MatrixType>& tri_) : tri(tri_) {}

    template<class MatrixType>
    template<class OtherMatrix>
    void TridiagonalMatrixT<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        const size_t order = getRow();
        target = RealType::Zero();
        target(0, 0) = tri.working.calc(0, 0);
        target(1, 0) = tri.normBuffer[0];
        size_t i = 1;
        for (; i < order - 2; ++i) {
            target(i - 1, i) = target(i, i - 1).conjugate();
            target(i, i) = tri.working.calc(i, i);
            target(i + 1, i) = tri.normBuffer[i];
        }
        target(i - 1, i) = target(i, i - 1).conjugate();
        target(i, i) = tri.working.calc(i, i);
        target(i + 1, i) = tri.working.calc(i + 1, i);
        ++i;
        target(i - 1, i) = target(i, i - 1).conjugate();
        target(i, i) = tri.working.calc(i, i);
    }
}
