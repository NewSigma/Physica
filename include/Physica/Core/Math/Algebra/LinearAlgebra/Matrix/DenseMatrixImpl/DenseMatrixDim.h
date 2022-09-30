/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    using Utils::Dynamic;

    template<class Derived, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim {
    public:
        DenseMatrixDim() = default;
        DenseMatrixDim([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) { assert(row_ == Row); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
        /* Operations */
        void resize([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) { /* Do nothing */ }
        /* Helper */
        void swap([[maybe_unused]] DenseMatrixDim& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim<Derived, Dynamic, Column, MaxRow, MaxColumn> : public Utils::CRTPBase<Derived, 2> {
        using Base = Utils::CRTPBase<Derived, 2>;
    public:
        DenseMatrixDim() = default;
        DenseMatrixDim([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) { /* Do nothing */ }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(size % getColumn() == 0);
            return size / getColumn();
        }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
        /* Operations */
        void resize([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) { /* Do nothing */ }
        /* Helper */
        void swap([[maybe_unused]] DenseMatrixDim& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim<Derived, Row, Dynamic, MaxRow, MaxColumn> : public Utils::CRTPBase<Derived, 2> {
        using Base = Utils::CRTPBase<Derived, 2>;
    public:
        DenseMatrixDim() = default;
        DenseMatrixDim([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) { assert(row_ == Row); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(size % getRow() == 0);
            return size / getRow();
        }
        /* Operations */
        void resize([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) { /* Do nothing */ }
        /* Helper */
        void swap([[maybe_unused]] DenseMatrixDim& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim<Derived, Dynamic, Dynamic, MaxRow, MaxColumn> : public Utils::CRTPBase<Derived, 2> {
        using Base = Utils::CRTPBase<Derived, 2>;
    private:
        size_t r;
    public:
        DenseMatrixDim() : r(0) {}
        DenseMatrixDim([[maybe_unused]] size_t row_, [[maybe_unused]] size_t column_) : r(row_) {}
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return r; }
        [[nodiscard]] size_t getColumn() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(r == 0 || size % getRow() == 0);
            return r == 0 ? 0 : size / getRow();
        }
        /* Operations */
        void resize(size_t row_, [[maybe_unused]] size_t column_) { r = row_; }
        /* Helper */
        void swap(DenseMatrixDim& dim) noexcept { std::swap(r, dim.r); }
    };
}
