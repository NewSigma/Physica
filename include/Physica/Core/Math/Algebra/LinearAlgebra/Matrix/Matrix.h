#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include <memory>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"

namespace Physica::Core {
    enum MatrixType {
        Row, Column
    };

    template<class T = MultiScalar, MatrixType type = Column, size_t maxRow = Dynamic, size_t maxColumn = Dynamic>
    class Matrix;
    /*!
     * Specialization for column matrix.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    class Matrix<T, Column, maxRow, maxColumn>
            : public CStyleArray<Vector<T, maxRow>, maxColumn> {
    public:
        typedef CStyleArray<Vector<T, maxRow>, maxColumn> Base;
        typedef Vector<T, maxRow> VectorType;

        Matrix() = default;
        explicit Matrix(size_t length);
        Matrix(const Matrix& matrix) = default;
        Matrix(Matrix&& matrix) noexcept;
        ~Matrix() = default;
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) { return Base::operator[](c)[r]; }
        [[nodiscard]] const T& operator()(size_t r, size_t c) const  { return Base::operator[](c)[r]; }
        Matrix& operator=(const Matrix& m) = default;
        Matrix& operator=(Matrix&& m) noexcept;
        /* Matrix Operations */
        void appendRow(const Vector<T, maxRow>& v);
        void appendRow(Vector<T, maxRow>&& v) noexcept;
        void appendColumn(const Vector<T, maxRow>& v);
        void appendColumn(Vector<T, maxRow>&& v) noexcept;
        void appendMatrixRow(const Matrix& m);
        void appendMatrixRow(Matrix&& m);
        void appendMatrixColumn(const Matrix& m);
        void appendMatrixColumn(Matrix&& m);
        void removeRowAt(size_t index);
        inline void removeColumnAt(size_t index);
        Vector<T, maxRow> cutRow();
        Vector<T, maxRow> cutColumn();
        Matrix<T, Column, Dynamic, Dynamic> cutMatrixRow(size_t from);
        Matrix<T, Column, Dynamic, Dynamic> cutMatrixColumn(size_t from);
        void rowSwap(size_t r1, size_t r2) noexcept;
        void columnSwap(size_t c1, size_t c2) noexcept;
        void rowReduce(size_t r1, size_t r2, size_t element);
        void columnReduce(size_t c1, size_t c2, size_t element);
        /* Getters */
        [[nodiscard]] constexpr static MatrixType getType() { return Column; }
        [[nodiscard]] size_t getRow() const { return Base::operator[](0).getLength(); }
        [[nodiscard]] size_t getColumn() const { return Base::getLength(); }
    };
    /*!
     * Specialization for row matrix.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    class Matrix<T, Row, maxRow, maxColumn>
            : public CStyleArray<Vector<T, maxColumn>, maxRow> {
    public:
        typedef CStyleArray<Vector<T, maxColumn>, maxRow> Base;
        typedef Vector<T, maxRow> VectorType;

        Matrix() = default;
        explicit Matrix(size_t length);
        Matrix(const Matrix& matrix) = default;
        Matrix(Matrix&& matrix) noexcept;
        ~Matrix() = default;
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) { return Base::operator[](r)[c]; }
        [[nodiscard]] const T& operator()(size_t r, size_t c) const  { return Base::operator[](r)[c]; }
        Matrix& operator=(const Matrix& m) = default;
        Matrix& operator=(Matrix&& m) noexcept;
        /* Matrix Operations */
        void appendRow(const Vector<T, maxColumn>& v);
        void appendRow(Vector<T, maxColumn>&& v) noexcept;
        void appendColumn(const Vector<T, maxColumn>& v);
        void appendColumn(Vector<T, maxColumn>&& v) noexcept;
        void appendMatrixRow(const Matrix& m);
        void appendMatrixRow(Matrix&& m);
        void appendMatrixColumn(const Matrix& m);
        void appendMatrixColumn(Matrix&& m);
        inline void removeRowAt(size_t index);
        void removeColumnAt(size_t index);
        Vector<T, maxColumn> cutRow();
        Vector<T, maxColumn> cutColumn();
        Matrix<T, Row, Dynamic, Dynamic> cutMatrixRow(size_t from);
        Matrix<T, Row, Dynamic, Dynamic> cutMatrixColumn(size_t from);
        void rowSwap(size_t r1, size_t r2) noexcept;
        void columnSwap(size_t c1, size_t c2) noexcept;
        void rowReduce(size_t r1, size_t r2, size_t element);
        void columnReduce(size_t c1, size_t c2, size_t element);
        /* Getters */
        [[nodiscard]] constexpr static MatrixType getType() { return Row; }
        [[nodiscard]] size_t getRow() const { return Base::getLength(); }
        [[nodiscard]] size_t getColumn() const { return Base::operator[](0).getLength(); }
    };
    /* Operators */
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    std::ostream& operator<<(std::ostream& os, const Matrix<T, type, maxRow, maxColumn>& m);

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> operator+(
            const Matrix<T, type, maxRow, maxColumn>& m1, const Matrix<T, type, maxRow, maxColumn>& m2);

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> operator-(
            const Matrix<T, type, maxRow, maxColumn>& m1, const Matrix<T, type, maxRow, maxColumn>& m2);

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> operator*(
            const Matrix<T, type, maxRow, maxColumn>& m1, const Matrix<T, type, maxRow, maxColumn>& m2);

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> operator*(const Matrix<T, type, maxRow, maxColumn>& m, const MultiScalar& n);
    /* Inline Implementations */
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    inline void operator+=(Matrix<T, type, maxRow, maxColumn>& m1
            , const Matrix<T, type, maxRow, maxColumn>& m2) { m1 = m1 + m2; }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    inline void operator-=(Matrix<T, type, maxRow, maxColumn>& m1
            , const Matrix<T, type, maxRow, maxColumn>& m2) { m1 = m1 - m2; }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    inline void operator*=(Matrix<T, type, maxRow, maxColumn>& m
            , const MultiScalar& n) { m = m * n; }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    inline void operator*=(Matrix<T, type, maxRow, maxColumn>& m1
            , const Matrix<T, type, maxRow, maxColumn>& m2) { m1 = m1 * m2; }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    inline void swap(Matrix<T, type, maxRow, maxColumn>& m1
            , Matrix<T, type, maxRow, maxColumn>& m2) noexcept { m1.swap(m2); }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> reciprocal(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> sqrt(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> factorial(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> ln(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> log(const Matrix<T, type, maxRow, maxColumn>& m, const MultiScalar& a);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> exp(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> pow(const Matrix<T, type, maxRow, maxColumn>& m, const MultiScalar& a);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> cos(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> sin(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> tan(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> sec(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> csc(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> cot(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arccos(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arcsin(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arctan(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arcsec(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arccsc(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arccot(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> cosh(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> sinh(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> tanh(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> sech(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> csch(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> coth(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arccosh(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arcsinh(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arctanh(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arcsech(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arccsch(const Matrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    Matrix<T, type, maxRow, maxColumn> arccoth(const Matrix<T, type, maxRow, maxColumn>& m);
}

#include "MatrixImpl.h"

#endif