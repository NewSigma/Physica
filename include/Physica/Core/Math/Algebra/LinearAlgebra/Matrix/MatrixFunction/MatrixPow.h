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
    template<Matrix T>
    class MatrixPow : public RValueMatrix<MatrixPow<T>> {
        using This = MatrixPow<T>;
        using Base = RValueMatrix<This>;

        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<T>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<T>();
        using TransposeRtnTy = std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using typename Base::ScalarType;
    private:
        const T& m;
        int power;
    public:
        MatrixPow(const T& m_, int power_);
        MatrixPow(const This&) = delete;
        MatrixPow(This&&) = delete;
        ~MatrixPow() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t, size_t) const { noImpl(); }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        /* Getters */
        [[nodiscard]] const T& getMatrix() const noexcept { return m; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return m.getCol(); }
        [[nodiscard]] int getPower() const noexcept { return power; }
    };

    template<Matrix T>
    MatrixPow<T>::MatrixPow(const T& m_, int power_) : m(m_), power(power_) {}

    template<Matrix T>
    [[nodiscard]] inline MatrixPow<T> pow(const T& m, int power) noexcept {
        return MatrixPow<T>(m, power);
    }
}

namespace Physica {
    template<Core::Matrix T>
    class Traits<Core::MatrixPow<T>> : public Traits<T> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
    };
}

#include "MatrixPowVecProduct.h"
