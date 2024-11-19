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
    template<Matrix T, Vector V>
    class MatrixVectorProduct;

    template<ExprType Type, Matrix T, class U, Vector V>
    class MatrixVectorProduct<MatrixExpr<Type, T, U>, V>
            : public RValueVector<MatrixVectorProduct<MatrixExpr<Type, T, U>, V>> {
        using MatrixType = MatrixExpr<Type, T, U>;
        using This = MatrixVectorProduct<MatrixType, V>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& expr;
        const V& vec;
    public:
        MatrixVectorProduct(const MatrixType& expr_, const V& vec_)
                : expr(expr_), vec(vec_) {
            assert(expr.getCol() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<LVector V1, class Executor = SequentialExecutor>
        inline void assignTo(V1& target) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return expr.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return expr; }
        [[nodiscard]] const V& getRHS() const noexcept { return vec; }
    private:
        template<LVector V1>
        inline void generalImpl(V1& target_) const;
    };

    template<ExprType Type, Matrix T, class U, Vector V>
    template<LVector V1, class Executor>
    inline void MatrixVectorProduct<MatrixExpr<Type, T, U>, V>::assignTo(V1& target) const {
        constexpr bool FastAssign = Traits<This>::FastAssign;
        if constexpr (!FastAssign)
            generalImpl(target);
        else if constexpr (Type == ExprType::Add) {
            using ExprAdd = decltype(expr.getLHS() * vec + expr.getRHS() * vec);
            target.template operator=<ExprAdd, Executor>(expr.getLHS() * vec + expr.getRHS() * vec);
        }
        else if constexpr (Type == ExprType::Sub) {
            using ExprSub = decltype(expr.getLHS() * vec - expr.getRHS() * vec);
            target.template operator=<ExprSub, Executor>(expr.getLHS() * vec - expr.getRHS() * vec);
        }
        else if constexpr (Type == ExprType::Mul) {
            using ExprMul = decltype((expr.getLHS() * vec) * expr.getRHS());
            target.template operator=<ExprMul, Executor>((expr.getLHS() * vec) * expr.getRHS());
        }
        else
            static_assert(!FastAssign, "[Error]: assignTo is not implemented");
    }

    template<ExprType Type, Matrix T, class U, Vector V>
    inline typename MatrixVectorProduct<MatrixExpr<Type, T, U>, V>::ScalarType
    MatrixVectorProduct<MatrixExpr<Type, T, U>, V>::calc(size_t index) const {
        return expr.row(index) * vec;
    }

    template<ExprType Type, Matrix T, class U, Vector V>
    template<LVector V1>
    inline void MatrixVectorProduct<MatrixExpr<Type, T, U>, V>::generalImpl(V1& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);

        if constexpr (is_continuous<V1>::value && Base::isReverseDiff)
            target.getDerived().makeContinuous();
    }
}

namespace Physica {
    template<ExprType Type, Matrix T, class U, Vector V>
    class Traits<MatrixVectorProduct<MatrixExpr<Type, T, U>, V>> {
        using MatrixType = MatrixExpr<Type, T, U>;
        using ExprType = ExprType;
        static_assert(MatrixType::ColAtCompile == V::SizeAtCompile || MatrixType::ColAtCompile == Dynamic || V::SizeAtCompile == Dynamic,
                      "Row and col do not match in matrix product");

        constexpr static bool calcFastAssign() {
            constexpr bool isScalarT2 = Scalar<U>;
            if constexpr (Type == ExprType::Add || Type == ExprType::Sub) {
                if constexpr (!isScalarT2) {
                    using LHS = decltype(std::declval<T>() * std::declval<V>());
                    using RHS = decltype(std::declval<U>() * std::declval<V>());
                    using ExprType = decltype(std::declval<LHS>() + std::declval<RHS>());
                    return Traits<ExprType>::FastAssign;
                }
                return false;
            }
            if constexpr (Type == ExprType::Mul) {
                return isScalarT2;
            }
            return false;
        }
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType, typename V::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;

        constexpr static bool FastAssign = calcFastAssign();
        constexpr static bool FastPacket = false;
    };
}
