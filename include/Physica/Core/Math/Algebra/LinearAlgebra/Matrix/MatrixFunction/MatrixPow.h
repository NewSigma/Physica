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
    template<class MatrixType> class MatrixPow;

    namespace Internal {
        template<class MatrixType>
        class Traits<MatrixPow<MatrixType>> : public Traits<MatrixType> {
        public:
            constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        };
    }

    template<class MatrixType>
    class MatrixPow : public RValueMatrix<MatrixPow<MatrixType>> {
        using This = MatrixPow<MatrixType>;
        using Base = RValueMatrix<This>;

        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixType& m;
        int power;
    public:
        MatrixPow(const RValueMatrix<MatrixType>& m_, int power_);
        MatrixPow(const This&) = delete;
        MatrixPow(This&&) = delete;
        ~MatrixPow() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t, size_t) const { throw NotImplementedException(); }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        /* Getters */
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return m; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return m.getColumn(); }
        [[nodiscard]] int getPower() const noexcept { return power; }
    };

    template<class MatrixType>
    MatrixPow<MatrixType>::MatrixPow(const RValueMatrix<MatrixType>& m_, int power_) : m(m_.getDerived()), power(power_) {}

    template<class MatrixType>
    [[nodiscard]] inline MatrixPow<MatrixType> pow(const RValueMatrix<MatrixType>& m, int power) noexcept {
        return MatrixPow<MatrixType>(m, power);
    }
}

#include "MatrixPowVecProduct.h"
