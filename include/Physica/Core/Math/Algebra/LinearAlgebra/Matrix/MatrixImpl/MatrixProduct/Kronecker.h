/*
 * Copyright 2025 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Scalar, size_t Order> class UnitMatrix;

    template<Matrix M1, Matrix M2>
    class Kronecker : public RValueMatrix<Kronecker<M1, M2>> {
        using This = Kronecker<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const M1& m1;
        const M2& m2;
    public:
        Kronecker(const M1& m1, const M2& m2);
        Kronecker(const This&) = default;
        Kronecker(This&&) noexcept = default;
        ~Kronecker() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_add(Matrix auto&& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
    };

    template<Matrix M1, Matrix M2>
    Kronecker<M1, M2>::Kronecker(const M1& m1, const M2& m2) : m1(m1), m2(m2) {}

    template<Matrix M1, Matrix M2>
    void Kronecker<M1, M2>::assign(Matrix auto&& target) const {
        for (size_t r = 0; r < m1.getRow(); ++r) {
            size_t offsetR = r * m2.getRow();
            if constexpr (instanceof_tx<UnitMatrix, M1>)
                m2.assign(target.block(offsetR, m2.getRow(), offsetR, m2.getCol()));
            else if constexpr (instanceof_tx<DiagMatrix, M2>)
                (m2 * m1.calc(r, r)).assign(target.block(offsetR, m2.getRow(), offsetR, m2.getCol()));
            else {
                for (size_t c = 0; c < m1.getCol(); ++c) {
                    size_t offsetC = c * m2.getCol();
                    (m2 * m1.calc(r, c)).assign(target.block(offsetR, m2.getRow(), offsetC, m2.getCol()));
                }
            }
        }
    }

    template<Matrix M1, Matrix M2>
    void Kronecker<M1, M2>::assign_add(Matrix auto&& target) const {
        for (size_t r = 0; r < m1.getRow(); ++r) {
            size_t offsetR = r * m2.getRow();
            if constexpr (instanceof_tx<UnitMatrix, M1>)
                m2.assign_add(target.block(offsetR, m2.getRow(), offsetR, m2.getCol()));
            else if constexpr (instanceof_tx<DiagMatrix, M2>)
                (m2 * m1.calc(r, r)).assign_add(target.block(offsetR, m2.getRow(), offsetR, m2.getCol()));
            else {
                for (size_t c = 0; c < m1.getCol(); ++c) {
                    size_t offsetC = c * m2.getCol();
                    (m2 * m1.calc(r, c)).assign_add(target.block(offsetR, m2.getRow(), offsetC, m2.getCol()));
                }
            }
        }
    }

    template<Matrix M1, Matrix M2>
    auto Kronecker<M1, M2>::calc(size_t row, size_t col) const -> T {
        size_t row1 = row / m2.getRow();
        size_t row2 = row % m2.getRow();
        size_t col1 = col / m2.getCol();
        size_t col2 = col % m2.getCol();
        return m1.calc(row1, col1) * m2.calc(row2, col2);
    }

    template<Matrix M1, Matrix M2>
    size_t Kronecker<M1, M2>::getRow() const noexcept {
        return m1.getRow() * m2.getRow();
    }
    
    template<Matrix M1, Matrix M2>
    size_t Kronecker<M1, M2>::getCol() const noexcept {
        return m1.getCol() * m2.getCol();
    }

    [[nodiscard]] auto kronecker(const Matrix auto& m1, const Matrix auto& m2) noexcept {
        return Kronecker(m1, m2);
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<Kronecker<M1, M2>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename M2::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor;
        constexpr static size_t RowAtCompile = M1::RowAtCompile * M2::RowAtCompile;
        constexpr static size_t ColAtCompile = M1::ColAtCompile * M2::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
