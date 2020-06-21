#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include <memory>
#include "Vector.h"

namespace Physica::Core {
    ////////////////////////////////////////Matrix////////////////////////////////////////////
    class Matrix : public CStyleArray<Vector, Dynamic> {
    public:
        enum MatrixType {
            Row, Column
        };
    protected:
        MatrixType type;
    public:
        explicit Matrix(MatrixType type);
        Matrix(size_t capacity, MatrixType type);
        Matrix(size_t length, size_t capacity, MatrixType type);
        explicit Matrix(const CStyleArray<Vector, Dynamic>& array, MatrixType type);
        explicit Matrix(CStyleArray<Vector, Dynamic>&& array, MatrixType type) noexcept;
        Matrix(const Matrix& matrix) = default;
        Matrix(Matrix&& matrix) noexcept;
        virtual ~Matrix() = default;
        /* Operators */
        [[nodiscard]] virtual Scalar& operator()(size_t row, size_t column) = 0;
        [[nodiscard]] virtual const Scalar& operator()(size_t row, size_t column) const = 0;
        Matrix& operator=(const Matrix& m);
        Matrix& operator=(Matrix&& m) noexcept;
        /* Matrix Operations */
        virtual void appendRow(const Vector& v) = 0;
        virtual void appendRow(Vector&& v) noexcept = 0;
        virtual void appendColumn(const Vector& v) = 0;
        virtual void appendColumn(Vector&& v) noexcept = 0;
        virtual void appendMatrixRow(const Matrix& m) = 0;
        virtual void appendMatrixRow(Matrix&& m) = 0;
        virtual void appendMatrixColumn(const Matrix& m) = 0;
        virtual void appendMatrixColumn(Matrix&& m) = 0;
        virtual Vector cutRow() = 0;
        virtual Vector cutColumn() = 0;
        virtual std::unique_ptr<Matrix> cutMatrixRow(size_t from) = 0;
        virtual std::unique_ptr<Matrix> cutMatrixColumn(size_t from) = 0;
        virtual void rowSwap(size_t r1, size_t r2) noexcept = 0;
        virtual void columnSwap(size_t c1, size_t c2) noexcept = 0;
        virtual void rowReduce(size_t r1, size_t r2, size_t element) = 0;
        virtual void columnReduce(size_t c1, size_t c2, size_t element) = 0;
        void upperEliminate(size_t index);
        void lowerEliminate(size_t index);
        void impactPivoting();
        /* Helpers */
        void swap(Matrix& m) noexcept;
        [[nodiscard]] MatrixType getType() const noexcept { return type; }
        [[nodiscard]] virtual size_t row() const = 0;
        [[nodiscard]] virtual size_t column() const = 0;
        /* Pivoting */
        inline size_t completePivoting(size_t column);
        inline size_t partialPivoting(size_t column);
        inline void impactPivoting(size_t row);
        /* Friends */
        friend class ColumnMatrix;
        friend class RowMatrix;
        friend class LinearEquations;
        friend class LUDecomposition;
    };
    /* Operators */
    std::ostream& operator<<(std::ostream& os, const Matrix& m);
    std::unique_ptr<Matrix> operator+(const Matrix& m1, const Matrix& m2);
    std::unique_ptr<Matrix> operator-(const Matrix& m1, const Matrix& m2);
    std::unique_ptr<Matrix> operator*(const Matrix& m, const Scalar& n);
    std::unique_ptr<Matrix> operator*(const Matrix& m1, const Matrix& m2);
    /* Inline Implementations */
    inline void operator+=(Matrix& m1, const Matrix& m2) { m1 = *(m1 + m2); }
    inline void operator-=(Matrix& m1, const Matrix& m2) { m1 = *(m1 - m2); }
    inline void operator*=(Matrix& m, const Scalar& n) { m = *(m * n); }
    inline void operator*=(Matrix& m1, const Matrix& m2) { m1 = *(m1 * m2); }
    /*!
     * Select the main element of a column of Matrix and execute a row swap as well as a column swap to
     * make it stands at (column, column), return the origin column index of the main element.
     * The return values should be stored to recover the correct order of the solution.
     *
     * Complexity: O((rank - column) ^ 2)
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.35
     */
    inline size_t Matrix::completePivoting(size_t column) {
        Q_ASSERT(column < this->column());
        const auto rank = row();
        size_t main_row_index = 0, main_column_index = 0;
        const Scalar* main = &BasicConst::getInstance().get_0();
        for(size_t i = column; i < rank; ++i) {
            for(size_t j = column; j < rank; ++j) {
                const Scalar* temp = &(*this)(i, j);
                bool larger = absCompare(*main, *temp);
                main = larger ? main : temp;
                main_row_index = larger ? main_row_index : j;
                main_column_index = larger ? main_column_index : i;
            }
        }
        this->rowSwap(column, main_row_index);
        return main_column_index;
    }
    /*!
     * Select the main element of a column of Matrix and execute row swaps to make its row index equals to column index,
     * return the origin row index of the main element.
     *
     * Complexity: O(rank - column)
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.35
     */
    inline size_t Matrix::partialPivoting(size_t column) {
        Q_ASSERT(column < this->column());
        const auto rank = row();
        size_t main_index = column;
        const Scalar* main = &(*this)(column, column);
        for(size_t j = column + 1; j < rank; ++j) {
            const Scalar* temp = &(*this)(j, column);
            bool larger = absCompare(*main, *temp);
            main = larger ? main : temp;
            main_index = larger ? main_index : j;
        }
        this->rowSwap(column, main_index);
        return main_index;
    }
    /*!
     * Divide the row by the element with the largest abstract value. \row is a row of a matrix.
     * \rank is the rank of the equations.
     *
     * Complexity: O(rank)
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.35
     */
    inline void Matrix::impactPivoting(size_t row) {
        Q_ASSERT(row < this->row());
        const auto col = column();
        const Scalar* main = &(*this)(row, 0);
        for(size_t i = 1; i < col; ++i) {
            const Scalar* temp = &(*this)(row, i);
            main = absCompare(*main, *temp) ? main : temp;
        }
        const Scalar n = reciprocal(*main);
        for(size_t i = 0; i < col; ++i)
            (*this)(row, i) *= n;
    }

    inline void swap(Matrix& m1, Matrix& m2) noexcept { m1.swap(m2); }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    std::unique_ptr<Matrix> reciprocal(const Matrix& n);
    std::unique_ptr<Matrix> sqrt(const Matrix& n);
    std::unique_ptr<Matrix> factorial(const Matrix& n);
    std::unique_ptr<Matrix> ln(const Matrix& n);
    std::unique_ptr<Matrix> log(const Matrix& n, const Scalar& a);
    std::unique_ptr<Matrix> exp(const Matrix& n);
    std::unique_ptr<Matrix> pow(const Matrix& n, const Scalar& a);
    std::unique_ptr<Matrix> cos(const Matrix& n);
    std::unique_ptr<Matrix> sin(const Matrix& n);
    std::unique_ptr<Matrix> tan(const Matrix& n);
    std::unique_ptr<Matrix> sec(const Matrix& n);
    std::unique_ptr<Matrix> csc(const Matrix& n);
    std::unique_ptr<Matrix> cot(const Matrix& n);
    std::unique_ptr<Matrix> arccos(const Matrix& n);
    std::unique_ptr<Matrix> arcsin(const Matrix& n);
    std::unique_ptr<Matrix> arctan(const Matrix& n);
    std::unique_ptr<Matrix> arcsec(const Matrix& n);
    std::unique_ptr<Matrix> arccsc(const Matrix& n);
    std::unique_ptr<Matrix> arccot(const Matrix& n);
    std::unique_ptr<Matrix> cosh(const Matrix& n);
    std::unique_ptr<Matrix> sinh(const Matrix& n);
    std::unique_ptr<Matrix> tanh(const Matrix& n);
    std::unique_ptr<Matrix> sech(const Matrix& n);
    std::unique_ptr<Matrix> csch(const Matrix& n);
    std::unique_ptr<Matrix> coth(const Matrix& n);
    std::unique_ptr<Matrix> arccosh(const Matrix& n);
    std::unique_ptr<Matrix> arcsinh(const Matrix& n);
    std::unique_ptr<Matrix> arctanh(const Matrix& n);
    std::unique_ptr<Matrix> arcsech(const Matrix& n);
    std::unique_ptr<Matrix> arccsch(const Matrix& n);
    std::unique_ptr<Matrix> arccoth(const Matrix& n);
}

#endif
