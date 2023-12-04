/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h"

namespace Physica::Core {
    /**
     * \class DenseMatrixExpression represents \param T1 \param type \param T2. e.g. matrix + scalar, expression * expression
     */
    template<ExpressionType type,
             class T1,
             class T2 = T1,
             class ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type>
    class DenseMatrixExpression;

    namespace Internal {
        template<ExpressionType type, class T1, class T2, class ResultType>
        class Traits<DenseMatrixExpression<type, T1, T2, ResultType>> {
            constexpr static bool SameMajor = MatrixOption::isSameMajor<T1, T2>();
            constexpr static int Major = SameMajor ? MatrixOption::getMajor<T1>()
                                                   : int(MatrixOption::AnyMajor);
            constexpr static bool SameStorage = MatrixOption::isSameStorage<T1, T2>();
            constexpr static int Storage = SameStorage ? MatrixOption::getStorage<T1>()
                                                       : int(MatrixOption::AnyStorage);
        public:
            using ScalarType = ResultType;
            constexpr static int Option = Major | Storage;
            //Optimize: T1 and T2 may not have same compiling size, for example, T1 may be fixed size and T2 may be dynamic
            constexpr static size_t RowAtCompile = T1::RowAtCompile;
            constexpr static size_t ColumnAtCompile = T1::ColumnAtCompile;
            constexpr static size_t MaxRowAtCompile = T1::MaxRowAtCompile;
            constexpr static size_t MaxColumnAtCompile = T1::MaxColumnAtCompile;
            constexpr static size_t SizeAtCompile = T1::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = T1::MaxSizeAtCompile;
        };

        template<ExpressionType type, class T1, class T2, class ResultType>
        class Traits<DenseMatrixExpression<type, T1, ScalarBase<T2>, ResultType>> {
        public:
            using ScalarType = ResultType;
            constexpr static int Option = T1::Option;
            constexpr static size_t RowAtCompile = T1::RowAtCompile;
            constexpr static size_t ColumnAtCompile = T1::ColumnAtCompile;
            constexpr static size_t MaxRowAtCompile = T1::MaxRowAtCompile;
            constexpr static size_t MaxColumnAtCompile = T1::MaxColumnAtCompile;
            constexpr static size_t SizeAtCompile = T1::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = T1::MaxSizeAtCompile;
        };

        template<class Derived>
        class DenseMatrixExpressionBase : public RValueMatrix<Derived> {
        public:
            using Base = RValueMatrix<Derived>;
            using typename Base::ScalarType;
        public:
            [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
            /* Getters */
            [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRow(); }
            [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
            [[nodiscard]] ScalarType sum() const {
                ScalarType result = 0;
                for (size_t major = 0; major < MatrixOption::selectMajor<Derived>(getRow(), getColumn()); ++major)
                    for (size_t minor = 0; minor < MatrixOption::selectMinor<Derived>(getRow(), getColumn()); ++minor)
                        result += Base::calcFromMajorMinor(major, minor);
                return result;
            }
        };
    }
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Minus, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Minus, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Minus, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_) : exp(exp_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return -exp.calc(row, col); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class DenseMatrixExpression<ExpressionType::Add, MatrixType1, MatrixType2>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Add, MatrixType1, MatrixType2>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Add, MatrixType1, MatrixType2>>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& exp1;
        const MatrixType2& exp2;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType1>& exp1_, const RValueMatrix<MatrixType2>& exp2_)
                : exp1(exp1_.getDerived()), exp2(exp2_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp1.calc(row, col)) + ScalarType(exp2.calc(row, col));
        }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp.calc(row, col)) + ScalarType(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class DenseMatrixExpression<ExpressionType::Sub, MatrixType1, MatrixType2>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sub, MatrixType1, MatrixType2>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sub, MatrixType1, MatrixType2>>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& exp1;
        const MatrixType2& exp2;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType1>& exp1_, const RValueMatrix<MatrixType2>& exp2_)
                : exp1(exp1_.getDerived()), exp2(exp2_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp1.calc(row, col)) - ScalarType(exp2.calc(row, col));
        }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(exp.calc(row, col)) - ScalarType(scalar); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class DenseMatrixExpression<ExpressionType::Mul, MatrixType1, MatrixType2>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Mul, MatrixType1, MatrixType2>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Mul, MatrixType1, MatrixType2>>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return mat1.calc(row, col) * mat2.calc(row, col);
        }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp.calc(row, col)) * ScalarType(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<ExpressionType::Div, MatrixType, AnyScalar>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Div, MatrixType, AnyScalar>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Div, MatrixType, AnyScalar>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(exp.calc(row, col)) / ScalarType(scalar); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Reciprocal,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Reciprocal, MatrixType, MatrixType, ScalarType>>;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return reciprocal(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Sqrt, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sqrt,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sqrt, MatrixType, MatrixType, ScalarType>>;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return sqrt(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Abs,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Abs, MatrixType, MatrixType, ScalarType>>;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return abs(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Square, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Square, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Square, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return square(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Ln, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Ln, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Ln, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ln(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Exp, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Exp, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Exp, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return exp(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Sin, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sin, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Sin, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return sin(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<ExpressionType::Cos, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Cos, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<ExpressionType::Cos, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return cos(mat.calc(row, col)); }
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    inline DenseMatrixExpression<ExpressionType::Minus, Derived>
    operator-(const RValueMatrix<Derived>& mat) {
        return DenseMatrixExpression<ExpressionType::Minus, Derived>(mat.getDerived());
    }

    template<ExpressionType type, class T1, class T2>
    inline DenseMatrixExpression<ExpressionType::Minus, DenseMatrixExpression<type, T1, T2>>
    operator-(const DenseMatrixExpression<type, T1, T2>& exp) {
        return DenseMatrixExpression<ExpressionType::Minus, DenseMatrixExpression<type, T1, T2>>(exp);
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>
    operator+(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>(mat, s);
    }

    template<class MatrixType1, class MatrixType2>
    inline DenseMatrixExpression<ExpressionType::Add, MatrixType1, MatrixType2>
    operator+(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        return DenseMatrixExpression<ExpressionType::Add, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<ExpressionType::Sub, MatrixType, ScalarBase<ScalarType>>
    operator-(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<ExpressionType::Sub, MatrixType, ScalarBase<ScalarType>>(mat, s);
    }

    template<class MatrixType1, class MatrixType2>
    inline DenseMatrixExpression<ExpressionType::Sub, MatrixType1, MatrixType2>
    operator-(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        return DenseMatrixExpression<ExpressionType::Sub, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const ScalarBase<ScalarType>& s, const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType1, class MatrixType2>
    inline DenseMatrixExpression<ExpressionType::Mul, MatrixType1, MatrixType2>
    hadamard(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        return DenseMatrixExpression<ExpressionType::Mul, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<ExpressionType::Div, MatrixType, ScalarType>
    operator/(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<ExpressionType::Div, MatrixType, ScalarType>(m, s);
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    reciprocal(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Sqrt, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    sqrt(const RValueMatrix<MatrixType>& m) {
        return {m};
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    abs(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Square, MatrixType> square(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Square, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Ln, MatrixType> ln(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Ln, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Exp, MatrixType> exp(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Exp, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Sin, MatrixType> sin(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Sin, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<ExpressionType::Cos, MatrixType> cos(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<ExpressionType::Cos, MatrixType>(m);
    }
}

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RealMatrix.h"
