/*
 * Copyright 2020-2024 Weibo He.
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

#include "../DenseSymmMatrix.h"

namespace Physica::Core {
    template<Scalar T, size_t Order, Vector U>
    class MatrixVectorProduct<DenseSymmMatrix<T, Order>, U> : public RValueVector<MatrixVectorProduct<DenseSymmMatrix<T, Order>, U>> {
        using MatrixType = DenseSymmMatrix<T, Order>;
        using This = MatrixVectorProduct<MatrixType, U>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    private:
        const MatrixType& mat;
        const U& vec;
    public:
        MatrixVectorProduct(const MatrixType& mat_, const U& vec_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& target) const;

        [[nodiscard]] inline CoDiff<ScalarType> calc(size_t index) const;
        [[nodiscard]] inline ValueType calc_value(size_t index) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Scalar T, size_t Order, Vector U>
    MatrixVectorProduct<DenseSymmMatrix<T, Order>, U>::MatrixVectorProduct(const MatrixType& mat_, const U& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Scalar T, size_t Order, Vector U>
    template<Vector V, class Executor>
    inline void MatrixVectorProduct<DenseSymmMatrix<T, Order>, U>::assign(V& target) const {
        const size_t length = getLength();
        assert(length == target.getLength());
        if (length >= 16) {
            target = T(0);
            for (size_t i = 0; i < length; ++i) {
                const size_t diag = mat.toIndex1D(i, i);
                const auto seg1 = mat.asVector().segment(diag, diag + length - i);
                const auto seg2 = vec.segment(i, length);
                target[i] += seg1 * seg2;

                if (i + 1 < length) {
                    auto seg = target.segment(i + 1, length);
                    seg += seg1.tail(1) * vec.calc(i);
                }
            }
        }
        else
            Base::assign(target);
    }

    template<Scalar T, size_t Order, Vector U>
    inline auto MatrixVectorProduct<DenseSymmMatrix<T, Order>, U>::calc(size_t index) const -> CoDiff<ScalarType> {
        return mat.row(index) * vec;
    }

    template<Scalar T, size_t Order, Vector U>
    inline auto MatrixVectorProduct<DenseSymmMatrix<T, Order>, U>::calc_value(size_t index) const -> ValueType {
        return mat.values().row(index) * vec.values();
    }
}

namespace Physica {
    template<Scalar T, size_t Order, Vector U>
    class Traits<MatrixVectorProduct<DenseSymmMatrix<T, Order>, U>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename U::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = U::SizeAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
