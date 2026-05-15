/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Scalar T>
    class HermiteToeplitzMatrix : public RValueMatrix<HermiteToeplitzMatrix<T>>{
        using This = HermiteToeplitzMatrix<T>;
        using Base = RValueMatrix<This>;
        // static_assert(T::isComplex(), "[Error]: Use SymmToeplitzMatrix instead"); // TODO: Add SymmToeplitzMatrix

        VectorND<T> firstCol;
    public:
        HermiteToeplitzMatrix() = default;
        explicit HermiteToeplitzMatrix(size_t size);
        explicit HermiteToeplitzMatrix(VectorND<T> firstCol_) noexcept;
        HermiteToeplitzMatrix(const This&) = default;
        HermiteToeplitzMatrix(This&&) noexcept = default;
        ~HermiteToeplitzMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const noexcept;

        [[nodiscard]] auto&& hermite(this auto&&) noexcept;
        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getFirstCol() const noexcept { return firstCol; }
        [[nodiscard]] size_t getRow() const noexcept { return firstCol.getLength(); }
        [[nodiscard]] size_t getCol() const noexcept { return firstCol.getLength(); }
    };

    template<Scalar T>
    HermiteToeplitzMatrix<T>::HermiteToeplitzMatrix(size_t size) : firstCol(size) {
        firstCol[0].imag() = 0;
    }

    template<Scalar T>
    HermiteToeplitzMatrix<T>::HermiteToeplitzMatrix(VectorND<T> firstCol_) noexcept : firstCol(std::move(firstCol_)) {
        assert(firstCol[0].imag().isZero() && "[Error]: Not a hermite matrix");
    }

    template<Scalar T>
    T HermiteToeplitzMatrix<T>::calc(size_t row, size_t col) const noexcept {
        bool flag = row >= col;
        size_t delta = flag ? row - col : col - row;
        T x = firstCol[delta];
        return flag ? x : x.conjugate();
    }

    template<Scalar T>
    auto&& HermiteToeplitzMatrix<T>::hermite(this auto&& self) noexcept {
        return std::forward<This>(self);
    }

    template<Scalar T>
    void HermiteToeplitzMatrix<T>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        firstCol.swap(obj.firstCol);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<HermiteToeplitzMatrix<T>> {
    public:
        using ScalarType = T;
        constexpr static int Major = MatrixMajor::BothMajor;
    };
}
