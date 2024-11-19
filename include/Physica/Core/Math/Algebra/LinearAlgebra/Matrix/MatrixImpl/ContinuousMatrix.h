/*
 * Copyright 2023 Weibo He.
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
    public:
        ~ContinuousMatrix() = default;
        /* Operators */
        using Base::operator=;
        inline This& operator=(const This& obj);
        This& operator=(This&& obj) noexcept = delete;
        Derived& operator=(const ScalarType& s);
        /* Operations */
        [[nodiscard]] inline RowVector row(size_t r);
        [[nodiscard]] inline const RowVector row(size_t r) const;
        [[nodiscard]] inline ColVector col(size_t c);
        [[nodiscard]] inline const ColVector col(size_t c) const;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline RowBlock<Row> rows(size_t fromRow, size_t rowCount);
        template<size_t Row = Dynamic>
        [[nodiscard]] inline const RowBlock<Row> rows(size_t fromRow, size_t rowCount) const;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline RowBlock<Row> topRows(size_t to);
        template<size_t Row = Dynamic>
        [[nodiscard]] inline const RowBlock<Row> topRows(size_t to) const;
        template<size_t Row = Dynamic>
        [[nodiscard]] inline RowBlock<Row> bottomRows(size_t from);
        template<size_t Row = Dynamic>
        [[nodiscard]] inline const RowBlock<Row> bottomRows(size_t from) const;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline ColBlock<Col> cols(size_t fromCol, size_t colCount);
        template<size_t Col = Dynamic>
        [[nodiscard]] inline const ColBlock<Col> cols(size_t fromCol, size_t colCount) const;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline ColBlock<Col> leftCols(size_t to);
        template<size_t Col = Dynamic>
        [[nodiscard]] inline const ColBlock<Col> leftCols(size_t to) const;
        template<size_t Col = Dynamic>
        [[nodiscard]] inline ColBlock<Col> rightCols(size_t from);
        template<size_t Col = Dynamic>
        [[nodiscard]] inline const ColBlock<Col> rightCols(size_t from) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> topLeftCorner(size_t toRow, size_t toCol);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> topLeftCorner(size_t toRow, size_t toCol) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> topLeftCorner(size_t to);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> topLeftCorner(size_t to) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> topRightCorner(size_t toRow, size_t fromCol);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> topRightCorner(size_t toRow, size_t fromCol) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> bottomLeftCorner(size_t fromRow, size_t toCol);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> bottomLeftCorner(size_t fromRow, size_t toCol) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> bottomRightCorner(size_t fromRow, size_t fromCol);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> bottomRightCorner(size_t fromRow, size_t fromCol) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> bottomRightCorner(size_t from);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> bottomRightCorner(size_t from) const;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline ContinuousMatrixBlock<Derived, Row, Col> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] inline const ContinuousMatrixBlock<Derived, Row, Col> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;

        void resize(size_t r, size_t c) { Base::getDerived().resize(r, c); }
        inline void makeContinuous();
        [[nodiscard]] bool checkContinuous() const;
        template<class RandomGenerator>
        void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);

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
    class Traits<Core::ContinuousMatrix<Derived>> : public Traits<Derived> {};
}

#include "ContinuousMatrixImpl/ContinuousMatrixImpl.h"
