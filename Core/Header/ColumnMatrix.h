/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COLUMNMATRIX_H
#define PHYSICA_COLUMNMATRIX_H

#include "Matrix.h"

namespace Physica::Core {
    class ColumnMatrix : virtual public Matrix {
    public:
        ColumnMatrix();
        ColumnMatrix(size_t length, size_t capacity);
        explicit ColumnMatrix(const CStyleArray<Vector>& array);
        explicit ColumnMatrix(CStyleArray<Vector>&& array) noexcept;
        ColumnMatrix(const ColumnMatrix& matrix) = default;
        ColumnMatrix(ColumnMatrix&& matrix) noexcept;
        /* Operators */
        [[nodiscard]] Numerical& operator()(size_t row, size_t column) override { return (*this)[column][row]; }
        [[nodiscard]] const Numerical& operator()(size_t row, size_t column) const override { return (*this)[column][row]; }
        /* Matrix Operations*/
        void appendRow(const Vector& v) override;
        void appendRow(Vector&& v) noexcept override;
        void appendColumn(const Vector& v) override;
        void appendColumn(Vector&& v) noexcept override;
        void appendMatrixRow(const Matrix& m) override;
        void appendMatrixRow(Matrix&& m) override;
        void appendMatrixColumn(const Matrix& m) override;
        void appendMatrixColumn(Matrix&& m) override;
        Vector cutRow() override;
        Vector cutColumn() override;
        std::unique_ptr<Matrix> cutMatrixRow(size_t from) override;
        std::unique_ptr<Matrix> cutMatrixColumn(size_t from) override;
        void rowSwap(size_t r1, size_t r2) noexcept override;
        void columnSwap(size_t c1, size_t c2) noexcept override;
        void rowReduce(size_t r1, size_t r2, size_t element) override;
        void columnReduce(size_t c1, size_t c2, size_t element) override;
        /* Helpers */
        [[nodiscard]] size_t row() const override { return (*this)[0].getLength(); }
        [[nodiscard]] size_t column() const override { return getLength(); }
    };
}

#endif