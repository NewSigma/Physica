/*
 * Copyright 2022-2024 Weibo He.
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

namespace Physica::Core {
    template<Scalar T>
    class PermutationMatrix : public RValueMatrix<PermutationMatrix<T>> {
        Array<size_t> indexes;
    public:
        PermutationMatrix() = default;
        PermutationMatrix(size_t order);
        PermutationMatrix(Array<size_t> indexes_);
        PermutationMatrix(const PermutationMatrix&) = default;
        PermutationMatrix(PermutationMatrix&&) noexcept = default;
        ~PermutationMatrix() = default;
        /* Operators */
        PermutationMatrix& operator=(PermutationMatrix obj) noexcept;
        /* Operations */
        void swapRows(size_t row1, size_t row2);
        [[nodiscard]] PermutationMatrix inverse() const noexcept;
        void swap(PermutationMatrix& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] size_t getRow() const noexcept { return indexes.getLength(); }
        [[nodiscard]] size_t getCol() const noexcept { return indexes.getLength(); }
    };

    template<Scalar T>
    PermutationMatrix<T>::PermutationMatrix(size_t order) : indexes(order) {
        for (size_t i = 0; i < order; ++i)
            indexes[i] = i;
    }

    template<Scalar T>
    PermutationMatrix<T>::PermutationMatrix(Array<size_t> indexes_) : indexes(std::move(indexes_)) {
        std::unordered_set<size_t> buffer{};
        for (size_t index : indexes) {
            if (index >= indexes.getLength())
                throw std::invalid_argument("[Error]: Invalid index");
            buffer.insert(index);
        }
        if (buffer.size() != indexes.getLength())
            throw std::invalid_argument("[Error]: Duplicate index is not allowed");
    }

    template<Scalar T>
    PermutationMatrix<T>& PermutationMatrix<T>::operator=(PermutationMatrix obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void PermutationMatrix<T>::swapRows(size_t row1, size_t row2) {
        std::swap(indexes[row1], indexes[row2]);
    }

    template<Scalar T>
    PermutationMatrix<T> PermutationMatrix<T>::inverse() const noexcept {
        Array<size_t> result(indexes.getLength());
        for (size_t i = 0; i < result.getLength(); ++i)
            result[indexes[i]] = i;
        return PermutationMatrix<T>(std::move(result));
    }

    template<Scalar T>
    void PermutationMatrix<T>::swap(PermutationMatrix& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        indexes.swap(obj.indexes);
    }

    template<Scalar T>
    T PermutationMatrix<T>::calc(size_t row, size_t col) const {
        return indexes[row] == col ? T(1) : T(0);
    }
}

namespace Physica {
    template<Core::Scalar T>
    class Traits<Core::PermutationMatrix<T>> {
    public:
        using ScalarType = T;
        constexpr static int Option = Core::MatrixOption::AnyMajor | Core::MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
