#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include "Vector.h"
#include <memory>

namespace Physica::Core {
    ////////////////////////////////////////Matrix////////////////////////////////////////////
    class Matrix {
    public:
        enum MatrixType {
            Row, Column
        };
    protected:
        Vector* __restrict vectors;
        size_t length;
        size_t capacity;
        MatrixType type;
    public:
        explicit Matrix(MatrixType type);
        Matrix(size_t length, MatrixType type);
        Matrix(const Matrix& matrix);
        Matrix(Matrix&& matrix) noexcept;
        virtual ~Matrix();
        /* Operators */
        [[nodiscard]] Vector& operator[](size_t index) { return vectors[index]; }
        [[nodiscard]] const Vector& operator[](size_t index) const { return vectors[index]; }
        [[nodiscard]] virtual Numerical& operator()(size_t row, size_t column) = 0;
        [[nodiscard]] virtual const Numerical& operator()(size_t row, size_t column) const = 0;
        Matrix& operator=(const Matrix& m) noexcept;
        Matrix& operator=(Matrix&& m) noexcept;
        /* Matrix Operations */
        virtual void appendRow(Vector v) = 0;
        virtual void appendColumn(Vector v) = 0;
        inline void appendMatrixRow(const Matrix& m);
        virtual void appendMatrixRow(Matrix&& m) = 0;
        inline void appendMatrixColumn(const Matrix& m);
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
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] size_t getCapacity() const noexcept { return capacity; }
        [[nodiscard]] MatrixType getType() const noexcept { return type; }
        [[nodiscard]] bool isEmpty() const noexcept { return length == 0; }
        [[nodiscard]] virtual size_t row() const = 0;
        [[nodiscard]] virtual size_t column() const = 0;
    protected:
        Matrix(Vector* vectors, size_t length, MatrixType type);
        void directSwap(size_t i, size_t j);
        void indirectSwap(size_t i, size_t j);
        inline void directReduce(size_t i, size_t j, size_t element);
        inline void indirectReduce(size_t i, size_t j, size_t element);
        inline void directVectorAppend(Vector&& v);
        inline void indirectVectorAppend(Vector&& v);
        inline void directMatrixAppend(Matrix&& m);
        inline void indirectMatrixAppend(Matrix&& m);
        inline Vector directVectorCut();
        inline Vector indirectVectorCut();
        inline Matrix* directMatrixCut(size_t from);
        inline Matrix* indirectMatrixCut(size_t from);
        /* Pivoting */
        inline size_t completePivoting(size_t column);
        inline void partialPivoting(size_t column);
        inline void impactPivoting(size_t row);
    private:
        /* Friends */
        friend class ColumnMatrix;
        friend class RowMatrix;
        friend class LinearEquations;
    };
    ////////////////////////////////////////Column Matrix////////////////////////////////////////////
    class ColumnMatrix : virtual public Matrix {
    public:
        ColumnMatrix();
        ColumnMatrix(size_t column, size_t row);
        ColumnMatrix(Vector* vectors, size_t length);
        /* Operators */
        [[nodiscard]] Numerical& operator()(size_t row, size_t column) override { return vectors[column][row]; }
        [[nodiscard]] const Numerical& operator()(size_t row, size_t column) const override { return vectors[column][row]; }
        /* Matrix Operations*/
        void appendRow(Vector v) override;
        void appendColumn(Vector v) override;
        void appendMatrixRow(Matrix&& m) override;
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
        [[nodiscard]] size_t row() const override { return vectors[0].getLength(); }
        [[nodiscard]] size_t column() const override { return length; }
    };
    ////////////////////////////////////////Row Matrix////////////////////////////////////////////
    class RowMatrix : virtual public Matrix {
    public:
        RowMatrix();
        RowMatrix(size_t column, size_t row);
        RowMatrix(Vector* vectors, size_t length);
        /* Operators */
        [[nodiscard]] Numerical& operator()(size_t row, size_t column) override { return vectors[row][column]; }
        [[nodiscard]] const Numerical& operator()(size_t row, size_t column) const override { return vectors[row][column]; }
        /* Matrix Operations */
        void appendRow(Vector v) override;
        void appendColumn(Vector v) override;
        void appendMatrixRow(Matrix&& m) override;
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
        [[nodiscard]] size_t row() const override { return length; }
        [[nodiscard]] size_t column() const override { return vectors[0].getLength(); }
    };
    /* Operators */
    std::ostream& operator<<(std::ostream& os, const Matrix& m);
    std::unique_ptr<Matrix> operator+(const Matrix& m1, const Matrix& m2);
    std::unique_ptr<Matrix> operator-(const Matrix& m1, const Matrix& m2);
    std::unique_ptr<Matrix> operator*(const Matrix& m, const Numerical& n);
    std::unique_ptr<Matrix> operator*(const Matrix& m1, const Matrix& m2);
    /* Inline Implementations */
    inline void operator+=(Matrix& m1, const Matrix& m2) { m1 = *(m1 + m2); }
    inline void operator-=(Matrix& m1, const Matrix& m2) { m1 = *(m1 - m2); }
    inline void operator*=(Matrix& m, const Numerical& n) { m = *(m * n); }
    inline void operator*=(Matrix& m1, const Matrix& m2) { m1 = *(m1 * m2); }

    inline void Matrix::appendMatrixRow(const Matrix& m) {
        appendMatrixRow(RowMatrix(reinterpret_cast<const RowMatrix&>(m))); //NOLINT assume m is a RowMatrix.
    }

    inline void Matrix::appendMatrixColumn(const Matrix& m) {
        appendMatrixColumn(ColumnMatrix(reinterpret_cast<const ColumnMatrix&>(m))); //NOLINT assume m is a RowMatrix.
    }

    inline void swap(Matrix& m1, Matrix& m2) noexcept { m1.swap(m2); }
    /*!
     * The following functions(direct and indirect) is the implementation of virtual function swaps and reduces.
     * Direct functions will not call Vector's member functions while indirect functions will.
     */
    inline void Matrix::directSwap(size_t i, size_t j) {
        Physica::Core::swap(vectors[i], vectors[j]);
    }

    inline void Matrix::indirectSwap(size_t i, size_t j) {
        for(size_t k = 0; k < length; ++k)
            Physica::Core::swap(vectors[i], vectors[j]);
    }
    //!Reduce the element at \j using \i.
    inline void Matrix::directReduce(size_t i, size_t j, size_t element) {
        auto& v1 = vectors[i];
        auto& v2 = vectors[j];
        v2 -= v1 * (v2[element] / v1[element]);
    }
    //!Reduce the element at \j using \i.
    inline void Matrix::indirectReduce(size_t i, size_t j, size_t element) {
        const auto dividend = vectors[element][i] / vectors[element][j];
        for(size_t k = 0; k < length; ++k)
            vectors[k][j] -= vectors[k][i] * dividend;
    }

    inline void Matrix::directVectorAppend(Vector&& v) {
        if(length == capacity) {
            ++capacity;
            vectors = reinterpret_cast<Vector*>(realloc(vectors, capacity * sizeof(Vector)));
        }
        new (vectors + length) Vector(std::move(v));
        ++length;
    }

    inline void Matrix::indirectVectorAppend(Vector&& v) {
        for(size_t i = 0; i < length; ++i)
            vectors[i].append(std::move(v[i]));
    }

    inline void Matrix::directMatrixAppend(Matrix&& m) {
        const auto new_length = length + m.getLength();
        if(new_length > capacity) {
            capacity = new_length;
            vectors = reinterpret_cast<Vector*>(realloc(vectors, capacity * sizeof(Vector)));
        }
        for(size_t i = length; i < new_length; ++i)
            new (vectors + i) Vector(std::move(m.vectors[i - length]));
        length = new_length;
    }

    inline void Matrix::indirectMatrixAppend(Matrix&& m) {
        for(size_t i = 0; i < length; ++i)
            vectors[i].append(std::move(m[i]));
    }

    inline Vector Matrix::directVectorCut() {
        Q_ASSERT(length > 0);
        --length;
        return Vector(std::move(vectors[length]));
    }

    inline Vector Matrix::indirectVectorCut() {
        Q_ASSERT(vectors[0].getLength() > 0);
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(vectors[i].cutLast());
        return result;
    }
    //!\from is included
    inline Matrix* Matrix::directMatrixCut(size_t from) {
        Q_ASSERT(from <= length);
        const auto new_length = length - from + 1;
        auto new_vectors = reinterpret_cast<Vector*>(malloc(new_length * sizeof(Vector)));
        length = from;
        for(size_t i = 0; i < new_length; ++i)
            new (new_vectors + i) Vector(std::move(vectors[from + i]));
        return getType() == Column ?
               static_cast<Matrix*>(new ColumnMatrix(new_vectors, new_length))
               : static_cast<Matrix*>(new RowMatrix(new_vectors, new_length));
    }
    //!\from is included
    inline Matrix* Matrix::indirectMatrixCut(size_t from) {
        Q_ASSERT(from <= vectors[0].getLength());
        auto new_vectors = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < length; ++i)
            new (new_vectors + i) Vector(vectors[i].cut(from));
        return getType() == Column ?
               static_cast<Matrix*>(new ColumnMatrix(new_vectors, length))
               : static_cast<Matrix*>(new RowMatrix(new_vectors, length));
    }
    /*!
     * Select the main element of a column of Matrix and execute a row swap as well as a column swap to
     * make it stands at (column, column), return the origin column index of the main element.
     * The return values should be stored to recover the correct order of the solution.
     *
     * Complexity: O((rank - column) ^ 2)
     */
    inline size_t Matrix::completePivoting(size_t column) {
        Q_ASSERT(column < this->column());
        const auto rank = row();
        auto& matrix = *this;
        size_t main_row_index = 0, main_column_index = 0;
        Numerical main(BasicConst::getInstance().get_0());
        for(size_t i = column; i < rank; ++i) {
            for(size_t j = column; j < rank; ++j) {
                Numerical temp = matrix(i, j);
                temp.toAbs();
                bool larger = main < temp;
                main = larger ? std::move(temp) : main;
                main_row_index = larger ? j : main_row_index;
                main_column_index = larger ? i : main_column_index;
            }
        }
        matrix.rowSwap(column, main_row_index);
        const Numerical reciprocal_main = reciprocal(main);
        for(size_t row = column; row < rank; ++row)
            matrix(row, column) *= reciprocal_main;
        return main_column_index;
    }
    /*!
     * Select the main element of a column of Matrix and execute row swaps to make its row index equals to column index,
     * return the origin row index of the main element.
     *
     * Complexity: O(rank - column)
     */
    inline void Matrix::partialPivoting(size_t column) {
        Q_ASSERT(column < this->column());
        const auto rank = row();
        auto& matrix = *this;
        size_t main_index = 0;
        Numerical main = matrix(column, column);
        main.toAbs();
        for(size_t j = column + 1; j < rank; ++j) {
            Numerical temp = matrix(j, column);
            temp.toAbs();
            bool larger = main < temp;
            main = larger ? std::move(temp) : main;
            main_index = larger ? j : main_index;
        }
        matrix.rowSwap(column, main_index);
        const Numerical reciprocal_main = reciprocal(main);
        for(size_t row = column; row < rank; ++row)
            matrix(row, column) *= reciprocal_main;
    }
    /*!
     * Divide the row by the element with the largest abstract value. Row is a row of a matrix.
     * Rank is the rank of the equations.
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
        auto& matrix = *this;
        Numerical main(matrix(row, 0));
        for(size_t i = 1; i < col; ++i) {
            Numerical temp(matrix(row, i));
            temp.toAbs();
            main = main < temp ? std::move(temp) : main;
        }
        main = reciprocal(main);
        for(size_t i = 0; i < col; ++i)
            matrix(row, i) *= main;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    std::unique_ptr<Matrix> reciprocal(const Matrix& n);
    std::unique_ptr<Matrix> sqrt(const Matrix& n);
    std::unique_ptr<Matrix> factorial(const Matrix& n);
    std::unique_ptr<Matrix> ln(const Matrix& n);
    std::unique_ptr<Matrix> log(const Matrix& n, const Numerical& a);
    std::unique_ptr<Matrix> exp(const Matrix& n);
    std::unique_ptr<Matrix> pow(const Matrix& n, const Numerical& a);
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
