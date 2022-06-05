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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h"

namespace Physica::Core {
    template<class MatrixType> class BiDiagMatrixB;

    namespace Internal {
        template<class MatrixType>
        class Traits<BiDiagMatrixB<MatrixType>> : public Traits<MatrixType> {
        private:
            using Base = Traits<MatrixType>;
            using Base::MatrixOption;
        };
    }
    /**
     * Decomposite matrix A like A = UBV^T
     * 
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.284-285
     */
    template<class MatrixType>
    class Bidiagonalization {
        using ScalarType = typename MatrixType::ScalarType;
        static_assert(!ScalarType::isComplex, "[Error]: Bidiagonalization do not support complex matrixes");
        using WorkingMatrix = typename MatrixType::ColMatrix;
        constexpr static size_t NumSingularValue = MatrixType::RowAtCompile > MatrixType::ColumnAtCompile
                                                                            ? MatrixType::ColumnAtCompile
                                                                            : MatrixType::RowAtCompile;
        using MainDiagVector = Vector<ScalarType, NumSingularValue>;
        using SubDiagVector = Vector<ScalarType, NumSingularValue == 0 ? Utils::Dynamic : NumSingularValue - 1>;
    private:
        WorkingMatrix working;
        MainDiagVector mainDiag;
        SubDiagVector subDiag;
    public:
        Bidiagonalization() = default;
        Bidiagonalization(size_t row, size_t col);
        Bidiagonalization(const MatrixType& source);
        Bidiagonalization(const Bidiagonalization&) = default;
        Bidiagonalization(Bidiagonalization&&) noexcept = default;
        ~Bidiagonalization() = default;
        /* Operators */
        Bidiagonalization& operator=(Bidiagonalization& obj) noexcept;
        /* Operations */
        template<class OtherMatrix>
        void compute(const RValueMatrix<OtherMatrix>& source);
        /* Getters */
        [[nodiscard]] BiDiagMatrixB<MatrixType> getMatrixB() const noexcept { return BiDiagMatrixB(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixU() const;
        [[nodiscard]] HouseholderSequence<Transpose<WorkingMatrix>> getMatrixV() const;
        /* Helpers */
        void swap(Bidiagonalization& obj) noexcept;
    private:
        void householderOnCol(size_t colIndex);
        friend class BiDiagMatrixB<MatrixType>;
    };

    template<class MatrixType>
    Bidiagonalization<MatrixType>::Bidiagonalization(size_t row, size_t col)
            : working(row, col)
            , mainDiag(col)
            , subDiag(col - 1) {
        assert(row >= col);
    }

    template<class MatrixType>
    Bidiagonalization<MatrixType>::Bidiagonalization(const MatrixType& source) : Bidiagonalization(source.getRow(), source.getColumn()) {
        compute(source);
    }

    template<class MatrixType>
    Bidiagonalization<MatrixType>& Bidiagonalization<MatrixType>::operator=(Bidiagonalization& obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class MatrixType>
    template<class OtherMatrix>
    void Bidiagonalization<MatrixType>::compute(const RValueMatrix<OtherMatrix>& source) {
        assert(source.getRow() >= source.getColumn());
        const size_t numCol = working.getColumn();
        /* Allocate */ {
            working = source;
            mainDiag.resize(numCol);
            subDiag.resize(numCol - 1);
        }

        size_t i = 0;
        ScalarType unit;
        for (; i < numCol - 2; ++i) {
            householderOnCol(i);

            auto row = working.row(i);
            auto sub_row = row.tail(i + 1);
            unit = sub_row[0].unit();
            subDiag[i] = -householderInPlace(sub_row) * unit;
            auto corner2 = working.bottomRightCorner(i + 1);
            applyHouseholder(corner2, sub_row);
        }
        /* Hangle last - 1 col */ {
            householderOnCol(i);
            subDiag[i] = working(i, i + 1);
            ++i;
        }
        /* Hangle last col */ {
            auto col = working.col(i);
            auto sub_col = col.tail(i);
            unit = sub_col[0].unit();
            mainDiag[i] = -householderInPlace(sub_col) * unit;
        }
    }

    template<class MatrixType>
    HouseholderSequence<typename Bidiagonalization<MatrixType>::WorkingMatrix>
    Bidiagonalization<MatrixType>::getMatrixU() const {
        HouseholderSequence result(working);
        result.setSize(working.getColumn());
        return result;
    }

    template<class MatrixType>
    HouseholderSequence<Transpose<typename Bidiagonalization<MatrixType>::WorkingMatrix>>
    Bidiagonalization<MatrixType>::getMatrixV() const {
        HouseholderSequence result(working.transpose());
        result.setSize(working.getColumn() - 2);
        result.setShift(1);
        return result;
    }

    template<class MatrixType>
    void Bidiagonalization<MatrixType>::swap(Bidiagonalization& obj) noexcept {
        working.swap(obj.working);
        mainDiag.swap(obj.mainDiag);
        subDiag.swap(obj.subDiag);
    }

    template<class MatrixType>
    void Bidiagonalization<MatrixType>::householderOnCol(size_t colIndex) {
        auto col = working.col(colIndex);
        auto sub_col = col.tail(colIndex);
        const ScalarType unit = sub_col[0].unit();
        mainDiag[colIndex] = -householderInPlace(sub_col) * unit;
        auto corner = working.bottomRightCorner(colIndex, colIndex + 1);
        applyHouseholder(sub_col, corner);
    }

    template<class MatrixType>
    inline void swap(Bidiagonalization<MatrixType>& obj1, Bidiagonalization<MatrixType>& obj2) noexcept {
        obj1.swap(obj2);
    }

    template<class MatrixType>
    class BiDiagMatrixB : public RValueMatrix<BiDiagMatrixB<MatrixType>> {
        using Base = RValueMatrix<BiDiagMatrixB<MatrixType>>;
        using typename Base::ScalarType;
        const Bidiagonalization<MatrixType>& bidiag;
    public:
        BiDiagMatrixB(const Bidiagonalization<MatrixType>& bidiag_) : bidiag(bidiag_) {}
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return bidiag.working.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return bidiag.working.getColumn(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void BiDiagMatrixB<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        target = ScalarType::Zero();
        const size_t col_1 = target.getColumn() - 1;
        for (size_t i = 0; i < col_1; ++i) {
            target(i, i) = bidiag.mainDiag[i];
            target(i, i + 1) = bidiag.subDiag[i];
        }
        target(col_1, col_1) = bidiag.mainDiag[col_1];
    }
}
