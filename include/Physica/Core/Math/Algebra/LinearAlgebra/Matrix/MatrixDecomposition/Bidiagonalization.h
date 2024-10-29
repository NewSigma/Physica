/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h>

namespace Physica::Core {
    template<class MatrixType> class BiDiagMatrixB;
    /**
     * Decomposite matrix A like A = UBV^T
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:284-285
     */
    template<class MatrixType>
    class Bidiagonalization {
        using This = Bidiagonalization<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
        using WorkingMatrix = typename MatrixType::ColMatrix;
        constexpr static size_t NumSingularValue = MatrixType::RowAtCompile > MatrixType::ColumnAtCompile
                                                                            ? MatrixType::ColumnAtCompile
                                                                            : MatrixType::RowAtCompile;
        using MainDiagVector = Vector<ScalarType, NumSingularValue>;
        using SubDiagVector = Vector<ScalarType, NumSingularValue == 0 ? Dynamic : NumSingularValue - 1>;

        static_assert(!ScalarType::isComplex, "[Error]: Bidiagonalization do not support complex matrixes");
    private:
        WorkingMatrix working;
        MainDiagVector mainDiag;
        SubDiagVector subDiag;
    public:
        Bidiagonalization() = default;
        Bidiagonalization(size_t row, size_t col);
        Bidiagonalization(const MatrixType& source);
        Bidiagonalization(const This&) = default;
        Bidiagonalization(This&&) noexcept = default;
        ~Bidiagonalization() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class OtherMatrix>
        void compute(const RValueMatrix<OtherMatrix>& source);
        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] BiDiagMatrixB<MatrixType> getMatrixB() const noexcept { return BiDiagMatrixB(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixU() const;
        [[nodiscard]] HouseholderSequence<WorkingMatrix, false> getMatrixV() const;
    private:
        void householderOnCol(size_t colIndex);
        friend class BiDiagMatrixB<MatrixType>;
    };

    template<class MatrixType>
    Bidiagonalization<MatrixType>::Bidiagonalization(size_t row, size_t col) {
        resize(row, col);
    }

    template<class MatrixType>
    Bidiagonalization<MatrixType>::Bidiagonalization(const MatrixType& source) : Bidiagonalization(source.getRow(), source.getColumn()) {
        compute(source);
    }

    template<class MatrixType>
    template<class OtherMatrix>
    void Bidiagonalization<MatrixType>::compute(const RValueMatrix<OtherMatrix>& source) {
        assert(source.getRow() >= source.getColumn());
        working = source;

        const size_t numCol = working.getColumn();
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
    void Bidiagonalization<MatrixType>::resize(size_t row, size_t col) {
        assert(row >= col && col > 1 && "[Error]: Invalid size");
        working.resize(row, col);
        mainDiag.resize(col);
        subDiag.resize(col - 1);
    }

    template<class MatrixType>
    void Bidiagonalization<MatrixType>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        mainDiag.swap(obj.mainDiag);
        subDiag.swap(obj.subDiag);
    }

    template<class MatrixType>
    HouseholderSequence<typename Bidiagonalization<MatrixType>::WorkingMatrix>
    Bidiagonalization<MatrixType>::getMatrixU() const {
        HouseholderSequence result(working);
        result.setSize(working.getColumn());
        return result;
    }

    template<class MatrixType>
    HouseholderSequence<typename Bidiagonalization<MatrixType>::WorkingMatrix, false>
    Bidiagonalization<MatrixType>::getMatrixV() const {
        HouseholderSequence<WorkingMatrix, false> result(working);
        result.setSize(working.getColumn() - 2);
        result.setShift(1);
        return result;
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
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return bidiag.working.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return bidiag.working.getColumn(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void BiDiagMatrixB<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        target = ScalarType(0);
        const size_t col_1 = target.getColumn() - 1;
        for (size_t i = 0; i < col_1; ++i) {
            target(i, i) = bidiag.mainDiag[i];
            target(i, i + 1) = bidiag.subDiag[i];
        }
        target(col_1, col_1) = bidiag.mainDiag[col_1];
    }
}

namespace Physica {
    template<class MatrixType>
    class Traits<Core::BiDiagMatrixB<MatrixType>> : public Traits<MatrixType> {
    private:
        using Base = Traits<MatrixType>;
        using Base::Option;
    };
}

namespace std {
    template<class MatrixType>
    inline void swap(Physica::Core::Bidiagonalization<MatrixType>& __restrict obj1, Physica::Core::Bidiagonalization<MatrixType>& __restrict obj2) noexcept {
        obj1.swap(obj2);
    }
}
