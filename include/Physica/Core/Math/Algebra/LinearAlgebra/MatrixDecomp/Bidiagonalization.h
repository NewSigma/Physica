/*
 * Copyright 2022-2025 Weibo He.
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

namespace Physica {
    template<Matrix M> class BiDiagMatrixB;
    /**
     * Decomposite matrix A like A = UBV^T
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:284-285
     */
    template<Matrix M>
    class Bidiagonalization {
        using This = Bidiagonalization<M>;
        using T = M::ScalarType;
        using WorkingMatrix = M::ColMatrix;
        constexpr static size_t NumSingularValue = M::RowAtCompile > M::ColAtCompile
                                                                            ? M::ColAtCompile
                                                                            : M::RowAtCompile;
        using MainDiagVector = DenseVector<T, NumSingularValue>;
        using SubDiagVector = DenseVector<T, NumSingularValue == 0 ? Dynamic : NumSingularValue - 1>;

        static_assert(!T::isComplex, "[Error]: Bidiagonalization do not support complex matrixes");
    private:
        WorkingMatrix working;
        MainDiagVector mainDiag;
        SubDiagVector subDiag;
    public:
        Bidiagonalization() = default;
        Bidiagonalization(size_t row, size_t col);
        Bidiagonalization(const M& source);
        Bidiagonalization(const This&) = default;
        Bidiagonalization(This&&) noexcept = default;
        ~Bidiagonalization() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] BiDiagMatrixB<M> getMatrixB() const noexcept { return BiDiagMatrixB(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixU() const;
        [[nodiscard]] HouseholderSequence<WorkingMatrix, false> getMatrixV() const;
    private:
        void householderOnCol(size_t colIndex);
        friend class BiDiagMatrixB<M>;
    };

    template<Matrix M>
    Bidiagonalization<M>::Bidiagonalization(size_t row, size_t col) {
        resize(row, col);
    }

    template<Matrix M>
    Bidiagonalization<M>::Bidiagonalization(const M& source) : Bidiagonalization(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<Matrix M>
    void Bidiagonalization<M>::compute(const Matrix auto& source) {
        assert(source.getRow() >= source.getCol());
        working = source;

        const size_t numCol = working.getCol();
        size_t i = 0;
        for (; i < numCol - 2; ++i) {
            householderOnCol(i);

            auto row = working.row(i);
            auto sub_row = row.tail(i + 1);
            auto factor = unit(sub_row[0]);
            subDiag[i] = -sub_row.householder() * factor;
            auto corner2 = working.bottomRightCorner(i + 1);
            applyHouseholder(corner2, sub_row);
        }
        // Handle last - 1 col
        householderOnCol(i);
        subDiag[i] = working(i, i + 1);
        ++i;
        // Handle last col
        if (working.getRow() != numCol) {
            auto col = working.col(i);
            auto sub_col = col.tail(i);
            auto factor = unit(sub_col[0]);
            mainDiag[i] = -sub_col.householder() * factor;
        }
        else {
            mainDiag[i] = -working(i, i);
            working(i, i) = 2;
        }
    }

    template<Matrix M>
    void Bidiagonalization<M>::resize(size_t row, size_t col) {
        assert(row >= col && col > 1 && "[Error]: Invalid size");
        working.resize(row, col);
        mainDiag.resize(col);
        subDiag.resize(col - 1);
    }

    template<Matrix M>
    void Bidiagonalization<M>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        mainDiag.swap(obj.mainDiag);
        subDiag.swap(obj.subDiag);
    }

    template<Matrix M>
    auto Bidiagonalization<M>::getMatrixU() const -> HouseholderSequence<WorkingMatrix> {
        HouseholderSequence result(working);
        result.setSize(working.getCol());
        return result;
    }

    template<Matrix M>
    auto Bidiagonalization<M>::getMatrixV() const -> HouseholderSequence<WorkingMatrix, false> {
        HouseholderSequence<WorkingMatrix, false> result(working);
        result.setSize(working.getCol() - 2);
        result.setShift(1);
        return result;
    }

    template<Matrix M>
    void Bidiagonalization<M>::householderOnCol(size_t colIndex) {
        auto col = working.col(colIndex);
        auto sub_col = col.tail(colIndex);
        const T factor = unit(sub_col[0]);
        mainDiag[colIndex] = -sub_col.householder() * factor;
        auto corner = working.bottomRightCorner(colIndex, colIndex + 1);
        applyHouseholder(sub_col, corner);
    }

    template<Matrix M>
    class BiDiagMatrixB : public RValueMatrix<BiDiagMatrixB<M>> {
        using Base = RValueMatrix<BiDiagMatrixB<M>>;
        using typename Base::T;

        const Bidiagonalization<M>& bidiag;
    public:
        BiDiagMatrixB(const Bidiagonalization<M>& bidiag_) : bidiag(bidiag_) {}
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return bidiag.working.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return bidiag.working.getCol(); }
    };

    template<Matrix M>
    void BiDiagMatrixB<M>::assign(Matrix auto& target) const {
        target = T(0);
        const size_t col_1 = target.getCol() - 1;
        for (size_t i = 0; i < col_1; ++i) {
            target(i, i) = bidiag.mainDiag[i];
            target(i, i + 1) = bidiag.subDiag[i];
        }
        target(col_1, col_1) = bidiag.mainDiag[col_1];
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<BiDiagMatrixB<M>> : public Traits<M> {
    private:
        using Base = Traits<M>;
        using Base::Option;
    };
}

namespace std {
    template<Physica::Matrix M>
    void swap(Physica::Bidiagonalization<M>& __restrict obj1, Physica::Bidiagonalization<M>& __restrict obj2) noexcept {
        obj1.swap(obj2);
    }
}
