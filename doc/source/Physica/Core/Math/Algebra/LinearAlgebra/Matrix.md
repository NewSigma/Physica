# Implementation details of matrix

Any matrix is divided into two parts:

MatrixImpl(contains data) -> MatrixBase(contains algorithm)

MatrixImpl extends MatrixBase so the algorithms are available to every matrix.

MatrixBase is a template and MatrixImpl is passed as template param so that the data is available to MatrixBase.

# API provided by Matrix

Each matrix should have a class Traits that declares the following public members:

ScalarType
int MatrixOption
size_t RowAtCompile
size_t ColumnAtCompile
size_t MaxRowAtCompile
size_t MaxColumnAtCompile
size_t SizeAtCompile
size_t MaxSizeAtCompile

MatrixBase uses APIs provided by MatrixImpl, the essential APIs are listed here:

> ScalarType& operator()(size_t row, size_t column);
> const ScalarType& operator()(size_t row, size_t column) const;
> (constexpr static) size_t getRow();
> (constexpr static) size_t getColumn();
