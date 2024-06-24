/*
 * Copyright 2024 WeiBo He.
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
    class MatrixVectorProduct<DenseMatrixExpression<Type, T1, T2, ResultType>, VectorType>
            : public RValueVector<MatrixVectorProduct<DenseMatrixExpression<Type, T1, T2, ResultType>, VectorType>> {
        using MatrixType = DenseMatrixExpression<Type, T1, T2, ResultType>;
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
        template<class OtherDerived>
        inline void assignTo(LValueVector<OtherDerived>& target_) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return expr.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return expr; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    template<class OtherDerived>
    inline void MatrixVectorProduct<DenseMatrixExpression<Type, T1, T2, ResultType>, VectorType>::assignTo(
            LValueVector<OtherDerived>& target_) const {
        constexpr bool isScalarT2 = is_scalar<T2>::value;
        auto& target = target_.getDerived();
        if constexpr (Type == ExpressionType::Mul && isScalarT2) {
            target = (expr.getMatrix() * vec) * expr.getScalar();
        }
        else {
            for (size_t i = 0; i < getLength(); ++i)
                target[i] = calc(i);
            
            constexpr bool isContinuous = std::is_base_of<ContinuousVector<OtherDerived>, OtherDerived>::value;
            if constexpr (isContinuous && Base::isReverseDiff)
                target.getDerived().makeContinuous();
        }
    }

    template<ExpressionType Type, class T1, class T2, class ResultType, class VectorType>
    inline typename MatrixVectorProduct<DenseMatrixExpression<Type, T1, T2, ResultType>, VectorType>::ScalarType
    MatrixVectorProduct<DenseMatrixExpression<Type, T1, T2, ResultType>, VectorType>::calc(size_t index) const {
        return expr.row(index) * vec;
    }
}
