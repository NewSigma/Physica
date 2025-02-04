/*
 * Copyright 2023-2025 Weibo He.
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

#include "LValueMatrix.h"
#include "ContinuousMatrixImpl/ContinuousMatrixBlock.h"
#include "ContinuousMatrixImpl/ContinuousFlatten.h"

namespace Physica::Core {
    /**
     * A ContinuousMatrix has its elements on major direction distributed continuously.
     */
    template<class Derived>
    class ContinuousMatrix : public LValueMatrix<Derived> {
        using Base = LValueMatrix<Derived>;
        using This = ContinuousMatrix<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
        using Base::isReverseDiff;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<This>();
        constexpr static bool isColMatrix = MatrixOption::isColMatrix<This>();
        using RowVector = std::conditional<isRowMatrix, ContinuousMatrixBlock<Derived, 1, ColAtCompile>, LMatrixBlock<Derived, 1, ColAtCompile>>::type;
        using ColVector = std::conditional<isColMatrix, ContinuousMatrixBlock<Derived, RowAtCompile, 1>, LMatrixBlock<Derived, RowAtCompile, 1>>::type;
        template<size_t Row>
        using RowBlock = ContinuousMatrixBlock<Derived, Row, ColAtCompile>;
        template<size_t Col>
        using ColBlock = ContinuousMatrixBlock<Derived, RowAtCompile, Col>;
        template<size_t Row, size_t Col>
        using BlockType = ContinuousMatrixBlock<Derived, Row, Col>;
    public:
        ~ContinuousMatrix() = default;
        /* Operators */
        using Base::operator=;
        inline This& operator=(const This& obj);
        This& operator=(This&& obj) noexcept = delete;
        Derived& operator=(const ScalarType& s);
        /* Operations */
        [[nodiscard]] auto toNumpy() const;

        [[nodiscard]] inline auto row(size_t r) noexcept;
        [[nodiscard]] inline const auto row(size_t r) const noexcept;
        [[nodiscard]] inline auto col(size_t c) noexcept;
        [[nodiscard]] inline const auto col(size_t c) const noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline auto rows(size_t fromRow, size_t rowCount) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline const auto rows(size_t fromRow, size_t rowCount) const noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline auto topRows(size_t to) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline const auto topRows(size_t to) const noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline auto bottomRows(size_t from) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline const auto bottomRows(size_t from) const noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline auto cols(size_t fromCol, size_t colCount) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline const auto cols(size_t fromCol, size_t colCount) const noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline auto leftCols(size_t to) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline const auto leftCols(size_t to) const noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline auto rightCols(size_t from) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline const auto rightCols(size_t from) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto topLeftCorner(size_t toRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto topLeftCorner(size_t toRow, size_t toCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto topLeftCorner(size_t to) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto topLeftCorner(size_t to) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto topRightCorner(size_t toRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto topRightCorner(size_t toRow, size_t fromCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto bottomLeftCorner(size_t fromRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto bottomRightCorner(size_t fromRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto bottomRightCorner(size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto bottomRightCorner(size_t from) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept;

        template<Matrix M>
        void resize(const M& m) { resize(m.getRow(), m.getCol()); }
        void resize(size_t r, size_t c) { Base::getDerived().resize(r, c); }

        template<RandomGenerator R>
        void random_uniform();
        template<RandomGenerator R>
        void random_normal();
        template<class Distribution, RandomGenerator R>
        void random_any(Distribution& dist);

        const H5DataSet<2> read(const H5Location& loc, const char* name);
        H5DataSet<2> write(H5Location& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ PtrTy data() { return Base::getDerived().data_ptr(0, 0); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data() const { return Base::getDerived().data_ptr(0, 0); }
        [[nodiscard]] ContinuousFlatten<Derived> flatten() { return {*this}; }
        [[nodiscard]] const ContinuousFlatten<Derived> flatten() const { return {const_cast<This&>(*this)}; }
    protected:
        ContinuousMatrix() = default;
        ContinuousMatrix(const This&) = default;
        ContinuousMatrix(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<ContinuousMatrix<Derived>> : public Traits<Derived> {};
}

#include "ContinuousMatrixImpl/ContinuousMatrixImpl.h" // IWYU pragma: export
