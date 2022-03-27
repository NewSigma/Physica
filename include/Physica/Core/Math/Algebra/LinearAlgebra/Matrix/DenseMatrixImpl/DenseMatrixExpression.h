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

#include "Physica/Utils/Template/ExpressionTemplateHelper.h"

namespace Physica::Core {
    /**
     * \class DenseMatrixExpression represents \param T1 \param type \param T2. e.g. matrix + scalar, expression * expression
     */
    template<Utils::ExpressionType type,
             class T1,
             class T2 = T1,
             class ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type>
    class DenseMatrixExpression;

    namespace Internal {
        template<Utils::ExpressionType type, class T1, class T2, class ResultType>
        class Traits<DenseMatrixExpression<type, T1, T2, ResultType>> {
            constexpr static bool SameMajor = DenseMatrixOption::isSameMajor<T1, T2>();
            constexpr static bool RowMajor = DenseMatrixOption::isRowMatrix<T1>();
            constexpr static int Major = SameMajor ? (RowMajor ? int(DenseMatrixOption::Column)
                                                               : int(DenseMatrixOption::Row))
                                                   : int(DenseMatrixOption::AnyMajor);
            constexpr static int Storage = (DenseMatrixOption::isElementMatrix<T1>() && DenseMatrixOption::isElementMatrix<T2>())
                                         ? DenseMatrixOption::Element
                                         : DenseMatrixOption::Vector;
        public:
            using ScalarType = ResultType;
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
            template<class OtherDerived>
            void assignTo(LValueMatrix<OtherDerived>& m) const {
                assert(getRow() == m.getRow() && getColumn() == m.getColumn());
                for (size_t i = 0; i < m.getMaxMajor(); ++i) {
                    for (size_t j = 0; j < m.getMaxMinor(); ++j) {
                        const size_t r = LValueMatrix<OtherDerived>::rowFromMajorMinor(i, j);
                        const size_t c = LValueMatrix<OtherDerived>::columnFromMajorMinor(i, j);
                        m(r, c) = calc(r, c);
                    }
                }       
            }

            [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
            /* Getters */
            [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getRow(); }
            [[nodiscard]] size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
        };
    }
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Minus, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Minus, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Minus, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_) : exp(exp_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return -exp.calc(row, col); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType1, MatrixType2>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType1, MatrixType2>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType1, MatrixType2>>;
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
        [[nodiscard]] size_t getRow() const { return exp1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>>;
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
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType1, MatrixType2>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType1, MatrixType2>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType1, MatrixType2>>;
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
        [[nodiscard]] size_t getRow() const { return exp1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(exp.calc(row, col)) - ScalarType(scalar); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType1, MatrixType2>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType1, MatrixType2>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType1, MatrixType2>>;
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
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>>;
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
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class MatrixType, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Div, MatrixType, AnyScalar>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Div, MatrixType, AnyScalar>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Div, MatrixType, AnyScalar>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(exp.calc(row, col)) / ScalarType(scalar); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Reciprocal,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Reciprocal, MatrixType, MatrixType, ScalarType>>;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return reciprocal(mat.calc(row, col)); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Abs,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Abs, MatrixType, MatrixType, ScalarType>>;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return abs(mat.calc(row, col)); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Square, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Square, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Square, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return square(mat.calc(row, col)); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Exp, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Exp, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Exp, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return exp(mat.calc(row, col)); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Sin, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Sin, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Sin, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return sin(mat.calc(row, col)); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class DenseMatrixExpression<Utils::ExpressionType::Cos, MatrixType>
            : public Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Cos, MatrixType>> {
    public:
        using Base = Internal::DenseMatrixExpressionBase<DenseMatrixExpression<Utils::ExpressionType::Cos, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        DenseMatrixExpression(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return cos(mat.calc(row, col)); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    inline DenseMatrixExpression<Utils::ExpressionType::Minus, Derived>
    operator-(const RValueMatrix<Derived>& mat) {
        return DenseMatrixExpression<Utils::ExpressionType::Minus, Derived>(mat.getDerived());
    }

    template<Utils::ExpressionType type, class T1, class T2>
    inline DenseMatrixExpression<Utils::ExpressionType::Minus, DenseMatrixExpression<type, T1, T2>>
    operator-(const DenseMatrixExpression<type, T1, T2>& exp) {
        return DenseMatrixExpression<Utils::ExpressionType::Minus, DenseMatrixExpression<type, T1, T2>>(exp);
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>
    operator+(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>(mat, s);
    }

    template<class MatrixType1, class MatrixType2>
    inline DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType1, MatrixType2>
    operator+(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        return DenseMatrixExpression<Utils::ExpressionType::Add, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType, ScalarBase<ScalarType>>
    operator-(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType, ScalarBase<ScalarType>>(mat, s);
    }

    template<class MatrixType1, class MatrixType2>
    inline DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType1, MatrixType2>
    operator-(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        return DenseMatrixExpression<Utils::ExpressionType::Sub, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const ScalarBase<ScalarType>& s, const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType1, class MatrixType2>
    inline DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType1, MatrixType2>
    hadamard(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        return DenseMatrixExpression<Utils::ExpressionType::Mul, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    inline DenseMatrixExpression<Utils::ExpressionType::Div, MatrixType, ScalarType>
    operator/(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) {
        return DenseMatrixExpression<Utils::ExpressionType::Div, MatrixType, ScalarType>(m, s);
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class MatrixType>
    DenseMatrixExpression<Utils::ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    reciprocal(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<Utils::ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    abs(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<Utils::ExpressionType::Square, MatrixType> square(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Square, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<Utils::ExpressionType::Exp, MatrixType> exp(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Exp, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<Utils::ExpressionType::Sin, MatrixType> sin(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Sin, MatrixType>(m);
    }

    template<class MatrixType>
    DenseMatrixExpression<Utils::ExpressionType::Cos, MatrixType> cos(const RValueMatrix<MatrixType>& m) {
        return DenseMatrixExpression<Utils::ExpressionType::Cos, MatrixType>(m);
    }
}
