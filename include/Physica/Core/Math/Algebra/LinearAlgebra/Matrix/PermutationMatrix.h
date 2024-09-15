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
    template<class ScalarType>
    class PermutationMatrix : public RValueMatrix<PermutationMatrix<ScalarType>> {
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
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] size_t getRow() const noexcept { return indexes.getLength(); }
        [[nodiscard]] size_t getColumn() const noexcept { return indexes.getLength(); }
    };

    template<class ScalarType>
    PermutationMatrix<ScalarType>::PermutationMatrix(size_t order) : indexes(order) {
        for (size_t i = 0; i < order; ++i)
            indexes[i] = i;
    }

    template<class ScalarType>
    PermutationMatrix<ScalarType>::PermutationMatrix(Array<size_t> indexes_) : indexes(std::move(indexes_)) {
        std::unordered_set<size_t> buffer{};
        for (size_t index : indexes) {
            if (index >= indexes.getLength())
                throw std::invalid_argument("[Error]: Invalid index");
            buffer.insert(index);
        }
        if (buffer.size() != indexes.getLength())
            throw std::invalid_argument("[Error]: Duplicate index is not allowed");
    }

    template<class ScalarType>
    PermutationMatrix<ScalarType>& PermutationMatrix<ScalarType>::operator=(PermutationMatrix obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    void PermutationMatrix<ScalarType>::swapRows(size_t row1, size_t row2) {
        std::swap(indexes[row1], indexes[row2]);
    }

    template<class ScalarType>
    PermutationMatrix<ScalarType> PermutationMatrix<ScalarType>::inverse() const noexcept {
        Array<size_t> result(indexes.getLength());
        for (size_t i = 0; i < result.getLength(); ++i)
            result[indexes[i]] = i;
        return PermutationMatrix<ScalarType>(std::move(result));
    }

    template<class ScalarType>
    void PermutationMatrix<ScalarType>::swap(PermutationMatrix& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        indexes.swap(obj.indexes);
    }

    template<class ScalarType>
    ScalarType PermutationMatrix<ScalarType>::calc(size_t row, size_t col) const {
        return indexes[row] == col ? ScalarType(1) : ScalarType(0);
    }
}

namespace Physica {
    using namespace Core;

    template<class T>
    class Traits<PermutationMatrix<T>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColumnAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
    };
}
