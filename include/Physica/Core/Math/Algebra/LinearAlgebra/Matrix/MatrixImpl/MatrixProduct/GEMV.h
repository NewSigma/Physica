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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    template<Matrix T, Vector U>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<T, U>> {
        using This = MatrixVectorProduct<T, U>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
    private:
        const T& mat;
        const U& vec;
    public:
        MatrixVectorProduct(const T& mat_, const U& vec_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& target) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const T& getLHS() const noexcept { return mat; }
        [[nodiscard]] const U& getRHS() const noexcept { return vec; }
    };

    template<Matrix T, Vector U>
    MatrixVectorProduct<T, U>::MatrixVectorProduct(const T& mat_, const U& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix T, Vector U>
    template<LVector V, class Executor>
    inline void MatrixVectorProduct<T, U>::assignTo(V& target) const {
        if constexpr (MatrixOption::isColMatrix<T>()) {
            target = mat.col(0) * vec.calc(0);
            for (size_t i = 1; i < vec.getLength(); ++i)
                target += mat.col(i) * vec.calc(i);
        }
        else {
            for (size_t i = 0; i < getLength(); ++i)
                target[i] = calc(i);
        }
    }

    template<Matrix T, Vector U>
    inline MatrixVectorProduct<T, U>::ScalarType MatrixVectorProduct<T, U>::calc(size_t index) const {
        return mat.row(index) * vec;
    }

    template<Matrix T, Vector U>
    [[nodiscard]] inline auto operator*(const T& mat, const U& vec) noexcept requires(T::RowAtCompile != 1) {
        return MatrixVectorProduct<T, U>(mat, vec);
    }

    template<Matrix T, Vector U>
    [[nodiscard]] inline auto operator*(const T& mat, const U& vec) requires(T::RowAtCompile == 1 && T::ColAtCompile == 1) {
        assert(mat.getCol() == vec.getLength());
        return mat.row(0) * vec;
    }
}

namespace Physica {
    template<Matrix T, Vector U>
    class Traits<Core::MatrixVectorProduct<T, U>> {
        static_assert(T::ColAtCompile == U::SizeAtCompile ||
                      T::ColAtCompile == Dynamic ||
                      U::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using ScalarType = Core::Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename U::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = T::RowAtCompile;
        constexpr static bool FastAssign = Core::MatrixOption::isColMatrix<T>();
        constexpr static bool FastPacket = false;
    };
}
