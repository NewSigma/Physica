/*
 * Copyright 2021-2022 WeiBo He.
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
    class HalfDenseMatrixBase {
        using ScalarType = typename Traits<Derived>::ScalarType;
        constexpr static size_t order = Traits<Derived>::RowAtCompile;
        constexpr static size_t maxOrder = Traits<Derived>::MaxRowAtCompile;
        constexpr static size_t size = order * (order + 1) / 2;
        constexpr static size_t maxSize = maxOrder * (maxOrder + 1) / 2;
    public:
        using Type = DenseMatrixStorageHelper<ScalarType, size, maxSize>;
    };

    template<class Derived>
    class HalfDenseMatrixStorage : public HalfDenseMatrixBase<Derived>::Type {
        constexpr static size_t MaxSizeAtCompile = Traits<Derived>::MaxSizeAtCompile;
        using T = typename Traits<Derived>::ScalarType;
        using Base = typename HalfDenseMatrixBase<Derived>::Type;
        size_t order;
    public:
        HalfDenseMatrixStorage() : Base(), order(MaxSizeAtCompile) {}
        HalfDenseMatrixStorage(size_t order_) : Base(order_ * (order_ + 1) / 2), order(order_) {}
        HalfDenseMatrixStorage(size_t order_, const T& t) : Base(order_ * (order_ + 1) / 2, t), order(order_) {}
        HalfDenseMatrixStorage(size_t row, [[maybe_unused]] size_t column) : HalfDenseMatrixStorage(row) { assert(row == column); }
        HalfDenseMatrixStorage(size_t row, size_t column, const T& t) : HalfDenseMatrixStorage(row, t) { assert(row == column); }
        HalfDenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
        /* Operators */
        using Base::operator[];
        /* Operations */
        void resize(size_t order_) { Base::resize(order_ * (order_ + 1) / 2); order = order_; }
        void resize(size_t row, [[maybe_unused]] size_t column) { assert(row == column); resize(row); order = row; }
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return order; }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getColumn() const noexcept { return getOrder(); }

        void swap(HalfDenseMatrixStorage& storage) noexcept;
    };

    template<class Derived>
    void HalfDenseMatrixStorage<Derived>::swap(HalfDenseMatrixStorage<Derived>& storage) noexcept {
        std::swap(order, storage.order);
        Base::swap(storage);
    }

    template<class Derived>
    inline void swap(Physica::Core::Internal::HalfDenseMatrixStorage<Derived>& mat1,
                     Physica::Core::Internal::HalfDenseMatrixStorage<Derived>& mat2) noexcept {
        mat1.swap(mat2);
    }
}
