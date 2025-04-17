/*
 * Copyright 2022-2025 Weibo He.
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

#include <unordered_set>
#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Scalar T>
    class PermutationMatrix : public RValueMatrix<PermutationMatrix<T>> {
        using This = PermutationMatrix<T>;

        Array<size_t> indices;
    public:
        PermutationMatrix() = default;
        PermutationMatrix(size_t order);
        PermutationMatrix(Array<size_t> indices_);
        PermutationMatrix(const This&) = default;
        PermutationMatrix(This&&) noexcept = default;
        ~PermutationMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swapRows(size_t row1, size_t row2);
        [[nodiscard]] PermutationMatrix inverse() const noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] auto& getIndices() noexcept { return indices; }
        [[nodiscard]] size_t getRow() const noexcept { return indices.getLength(); }
        [[nodiscard]] size_t getCol() const noexcept { return indices.getLength(); }
    };

    template<Scalar T>
    PermutationMatrix<T>::PermutationMatrix(size_t order) : indices(order) {
        for (size_t i = 0; i < order; ++i)
            indices[i] = i;
    }

    template<Scalar T>
    PermutationMatrix<T>::PermutationMatrix(Array<size_t> indices_) : indices(std::move(indices_)) {
        std::unordered_set<size_t> buffer{};
        for (size_t index : indices) {
            if (index >= indices.getLength())
                throw std::invalid_argument("[Error]: Invalid index");
            buffer.insert(index);
        }
        if (buffer.size() != indices.getLength())
            throw std::invalid_argument("[Error]: Duplicate index is not allowed");
    }

    template<Scalar T>
    void PermutationMatrix<T>::swapRows(size_t row1, size_t row2) {
        std::swap(indices[row1], indices[row2]);
    }

    template<Scalar T>
    PermutationMatrix<T> PermutationMatrix<T>::inverse() const noexcept {
        Array<size_t> result(indices.getLength());
        for (size_t i = 0; i < result.getLength(); ++i)
            result[indices[i]] = i;
        return PermutationMatrix<T>(std::move(result));
    }

    template<Scalar T>
    void PermutationMatrix<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        indices.swap(obj.indices);
    }

    template<Scalar T>
    T PermutationMatrix<T>::calc(size_t row, size_t col) const {
        return indices[row] == col ? T(1) : T(0);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<PermutationMatrix<T>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
