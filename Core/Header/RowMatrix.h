/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ROWMATRIX_H
#define PHYSICA_ROWMATRIX_H

#include "Matrix.h"

namespace Physica::Core {
    class RowMatrix : virtual public Matrix {
    public:
        RowMatrix();
        explicit RowMatrix(size_t capacity);
        RowMatrix(size_t length, size_t capacity);
        explicit RowMatrix(const CStyleArray<Vector, Dynamic>& array);
        explicit RowMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept;
        RowMatrix(const RowMatrix& matrix) = default;
        RowMatrix(RowMatrix&& matrix) noexcept;
        /* Operators */
        [[nodiscard]] Numerical& operator()(size_t row, size_t column) override { return (*this)[row][column]; }
        [[nodiscard]] const Numerical& operator()(size_t row, size_t column) const override { return (*this)[row][column]; }
        /* Matrix Operations */
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
        [[nodiscard]] size_t row() const override { return getLength(); }
        [[nodiscard]] size_t column() const override { return (*this)[0].getLength(); }
    };
}

#endif