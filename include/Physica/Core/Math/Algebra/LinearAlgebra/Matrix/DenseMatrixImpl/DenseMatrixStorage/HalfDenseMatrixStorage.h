/*
 * Copyright 2021 WeiBo He.
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

#include "AbstractDenseMatrixStorage.h"

namespace Physica::Core::Internal {
    template<class Derived>
    class HalfDenseMatrixStorage : public DenseMatrixStorageHelper<typename Traits<Derived>::ScalarType,
                                                                   Traits<Derived>::SizeAtCompile,
                                                                   Traits<Derived>::MaxSizeAtCompile> {
        using T = typename Traits<Derived>::ScalarType;
        using Base = DenseMatrixStorageHelper<T,
                                              Traits<Derived>::SizeAtCompile,
                                              Traits<Derived>::MaxSizeAtCompile>;
        size_t order;
    public:
        HalfDenseMatrixStorage() = default;
        HalfDenseMatrixStorage(size_t order_) : Base(order_ * (order_ + 1) / 2), order(order_) {}
        HalfDenseMatrixStorage(size_t order_, const T& t) : Base(order_ * (order_ + 1) / 2, t), order(order_) {}
        HalfDenseMatrixStorage(size_t row, [[maybe_unused]] size_t column) : HalfDenseMatrixStorage(row) { assert(row == column); }
        HalfDenseMatrixStorage(size_t row, size_t column, const T& t) : HalfDenseMatrixStorage(row, t) { assert(row == column); }
        HalfDenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
        /* Operators */
        using Base::operator[];
        /* Operations */
        void resize(size_t order_) { Base::resize(order_ * (order_ + 1) / 2); order = order_; }
        void resize(size_t row, [[maybe_unused]] size_t column) { assert(row == column); resize(row); }
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return order; }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getColumn() const noexcept { return getOrder(); }

        void swap(HalfDenseMatrixStorage& storage) { Base::swap(storage); }
    };
}
