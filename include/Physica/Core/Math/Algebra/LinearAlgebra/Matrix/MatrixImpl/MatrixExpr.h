/*
 * Copyright 2021-2024 WeiBo He.
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

#include <Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h>

namespace Physica::Core {
    /**
     * \class MatrixExpr represents \param T1 \param type \param T2. e.g. matrix + scalar, expression * expression
     */
    template<ExpressionType Type,
             class T1,
             class T2 = T1,
             class ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type>
    class MatrixExpr;
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class MatrixType>
    class MatrixExpr<ExpressionType::Minus, MatrixType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Minus, MatrixType>> {
        using This = MatrixExpr<ExpressionType::Minus, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& exp_) : exp(exp_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return -exp.calc(row, col); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>
            : public RValueMatrix<MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>> {
        using This = MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType1>() && MatrixOption::isSymmMatrix<MatrixType2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType1>() && MatrixOption::isHermiteMatrix<MatrixType2>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixExpr(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(mat1.calc(row, col)) + ScalarType(mat2.calc(row, col));
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat2.getColumn(); }
        [[nodiscard]] const MatrixType1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const MatrixType2& getRHS() const noexcept { return mat2; }
    };

    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>
            : public RValueMatrix<MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>> {
        using This = MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp.calc(row, col)) + ScalarType(scalar);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class MatrixExpr<ExpressionType::Sub, MatrixType1, MatrixType2>
            : public RValueMatrix<MatrixExpr<ExpressionType::Sub, MatrixType1, MatrixType2>> {
        using This = MatrixExpr<ExpressionType::Sub, MatrixType1, MatrixType2>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType1>() && MatrixOption::isSymmMatrix<MatrixType2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType1>() && MatrixOption::isHermiteMatrix<MatrixType2>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixExpr(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(mat1.calc(row, col)) - ScalarType(mat2.calc(row, col));
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat1.getColumn(); }
        [[nodiscard]] const MatrixType1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const MatrixType2& getRHS() const noexcept { return mat2; }
    };

    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>
            : public RValueMatrix<MatrixExpr<ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>> {
        using This = MatrixExpr<ExpressionType::Sub, MatrixType, ScalarBase<AnyScalar>>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(exp.calc(row, col)) - ScalarType(scalar); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    class MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>
            : public RValueMatrix<MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>> {
        using This = MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixExpr(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return mat1.calc(row, col) * mat2.calc(row, col);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>
            : public RValueMatrix<MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>> {
        using This = MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp.calc(row, col)) * ScalarType(scalar);
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return exp; }
        [[nodiscard]] const AnyScalar& getScalar() const noexcept { return scalar; }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExpressionType::Div, MatrixType, AnyScalar>
            : public RValueMatrix<MatrixExpr<ExpressionType::Div, MatrixType, AnyScalar>> {
        using This = MatrixExpr<ExpressionType::Div, MatrixType, AnyScalar>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(exp.calc(row, col)) / ScalarType(scalar); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class MatrixType>
    class MatrixExpr<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Reciprocal,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
        using This = MatrixExpr<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return reciprocal(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Sqrt, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Sqrt,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
        using This = MatrixExpr<ExpressionType::Sqrt, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return sqrt(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Abs,
                                                                               MatrixType,
                                                                               MatrixType,
                                                                               typename MatrixType::ScalarType::RealType>> {
        using This = MatrixExpr<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return abs(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Square, MatrixType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Square, MatrixType>> {
        using This = MatrixExpr<ExpressionType::Square, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return square(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Ln, MatrixType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Ln, MatrixType>> {
        using This = MatrixExpr<ExpressionType::Ln, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ln(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Exp, MatrixType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Exp, MatrixType>> {
        using This = MatrixExpr<ExpressionType::Exp, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return exp(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Sin, MatrixType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Sin, MatrixType>> {
        using This = MatrixExpr<ExpressionType::Sin, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return sin(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class MatrixExpr<ExpressionType::Cos, MatrixType>
            : public RValueMatrix<MatrixExpr<ExpressionType::Cos, MatrixType>> {
        using This = MatrixExpr<ExpressionType::Cos, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return cos(mat.calc(row, col)); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Minus, Derived>
    operator-(const RValueMatrix<Derived>& mat) noexcept {
        return MatrixExpr<ExpressionType::Minus, Derived>(mat.getDerived());
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>
    operator+(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        return {mat1, mat2};
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>
    operator+(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) noexcept {
        return {mat, s};
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Sub, MatrixType, ScalarBase<ScalarType>>
    operator-(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) noexcept {
        return MatrixExpr<ExpressionType::Sub, MatrixType, ScalarBase<ScalarType>>(mat, s);
    }

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Sub, MatrixType1, MatrixType2>
    operator-(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        return MatrixExpr<ExpressionType::Sub, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const ScalarBase<ScalarType>& s, const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) noexcept {
        return MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>
    hadamard(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        return MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>(mat1, mat2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Div, MatrixType, ScalarType>
    operator/(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) noexcept {
        return MatrixExpr<ExpressionType::Div, MatrixType, ScalarType>(m, s);
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    reciprocal_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Reciprocal, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>(m);
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Sqrt, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    sqrt_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return {m};
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>
    abs_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Abs, MatrixType, MatrixType, typename MatrixType::ScalarType::RealType>(m);
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Square, MatrixType> square_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Square, MatrixType>(m);
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Ln, MatrixType> ln_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Ln, MatrixType>(m);
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Exp, MatrixType> exp_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Exp, MatrixType>(m);
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Sin, MatrixType> sin_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Sin, MatrixType>(m);
    }

    template<class MatrixType>
    [[nodiscard]] MatrixExpr<ExpressionType::Cos, MatrixType> cos_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Cos, MatrixType>(m);
    }
}

namespace Physica {
    using namespace Core;

    template<ExpressionType Type, class T1, class T2, class ResultType>
    class Traits<MatrixExpr<Type, T1, T2, ResultType>> {
        constexpr static bool SameMajor = MatrixOption::isSameMajor<T1, T2>();
        constexpr static int Major = SameMajor ? MatrixOption::getMajor<T1>()
                                               : int(MatrixOption::AnyMajor);
        constexpr static bool SameStorage = MatrixOption::isSameStorage<T1, T2>();
        constexpr static int Storage = SameStorage ? MatrixOption::getStorage<T1>()
                                                   : int(MatrixOption::AnyStorage);
    public:
        using ScalarType = ResultType;
        constexpr static int Option = Major | Storage;
        // Optimize: T1 and T2 may not have same compiling size, for example, T1 may be fixed size and T2 may be dynamic
        constexpr static size_t RowAtCompile = T1::RowAtCompile;
        constexpr static size_t ColumnAtCompile = T1::ColumnAtCompile;
        constexpr static size_t MaxRowAtCompile = T1::MaxRowAtCompile;
        constexpr static size_t MaxColumnAtCompile = T1::MaxColumnAtCompile;
        constexpr static size_t SizeAtCompile = T1::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = T1::MaxSizeAtCompile;
    };

    template<ExpressionType Type, class T1, class T2, class ResultType>
    class Traits<MatrixExpr<Type, T1, ScalarBase<T2>, ResultType>> {
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
}

#include "MatrixExprImpl/ExprVecProduct.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixConvert.h"
