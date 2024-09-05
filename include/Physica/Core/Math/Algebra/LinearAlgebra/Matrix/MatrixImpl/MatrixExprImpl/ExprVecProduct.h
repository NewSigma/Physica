/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    template<class MatrixType, class VectorType> class MatrixVectorProduct;

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    class MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>
            : public RValueVector<MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>> {
        using MatrixType = MatrixExpr<Type, T1, T2, ResultType>;
        using This = MatrixVectorProduct<MatrixType, VectorType>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& expr;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const MatrixType& expr_, const RValueVector<VectorType>& vec_)
                : expr(expr_), vec(vec_.getDerived()) {
            assert(expr.getColumn() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& target_) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return expr.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return expr; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    private:
        template<class OtherDerived>
        inline void generalImpl(OtherDerived& target_) const;
    };

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    template<class OtherDerived, class Executor>
    inline void MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>::assignTo(
            LValueVector<OtherDerived>& target_) const {
        constexpr bool FastAssign = Traits<This>::FastAssign;
        auto& target = target_.getDerived();
        if constexpr (!FastAssign)
            generalImpl(target);
        else if constexpr (Type == ExpressionType::Add) {
            using ExprType = decltype(expr.getLHS() * vec + expr.getRHS() * vec);
            target.template operator=<ExprType, Executor>(expr.getLHS() * vec + expr.getRHS() * vec);
        }
        else if constexpr (Type == ExpressionType::Sub) {
            using ExprType = decltype(expr.getLHS() * vec - expr.getRHS() * vec);
            target.template operator=<ExprType, Executor>(expr.getLHS() * vec - expr.getRHS() * vec);
        }
        else if constexpr (Type == ExpressionType::Mul) {
            using ExprType = decltype((expr.getMatrix() * vec) * expr.getScalar());
            target.template operator=<ExprType, Executor>((expr.getMatrix() * vec) * expr.getScalar());
        }
        else
            static_assert(!FastAssign, "[Error]: assignTo is not implemented");

    }

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    inline typename MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>::ScalarType
    MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>::calc(size_t index) const {
        return expr.row(index) * vec;
    }

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    template<class OtherDerived>
    inline void MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>::generalImpl(OtherDerived& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);

        constexpr bool isContinuous = std::is_base_of<ContinuousVector<OtherDerived>, OtherDerived>::value;
        if constexpr (isContinuous && Base::isReverseDiff)
            target.getDerived().makeContinuous();
    }
}

namespace Physica {
    using namespace Core;

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    class Traits<MatrixVectorProduct<MatrixExpr<Type, T1, T2, ResultType>, VectorType>> {
        using MatrixType = MatrixExpr<Type, T1, T2, ResultType>;
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColumnAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");

        constexpr static bool calcFastAssign() {
            constexpr bool isScalarT2 = is_scalar<T2>::value;
            if constexpr (Type == ExpressionType::Add || Type == ExpressionType::Sub) {
                if constexpr (!isScalarT2) {
                    using LHS = decltype(std::declval<T1>() * std::declval<VectorType>());
                    using RHS = decltype(std::declval<T2>() * std::declval<VectorType>());
                    using ExprType = decltype(std::declval<LHS>() + std::declval<RHS>());
                    return Traits<ExprType>::FastAssign;
                }
                return false;
            }
            if constexpr (Type == ExpressionType::Mul) {
                return isScalarT2;
            }
            return false;
        }
    public:
        using ScalarType = ResultType;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;

        constexpr static bool FastAssign = calcFastAssign();
        constexpr static bool FastPacket = false;
    };
}
