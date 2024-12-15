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

#include "../DiffDenseMatrix.h"

namespace Physica::Core {
#define tparams Scalar T, DiffMode Mode, int Order, int Option, size_t Row, size_t Col
#define DiffDenseMatrix DenseMatrix<Diff<T, Mode, Order>, Option, Row, Col>

    template<tparams>
    DiffDenseMatrix::DenseMatrix(size_t row, size_t col) : v(row, col), g(row, col) {}

    template<tparams>
    DiffDenseMatrix::DenseMatrix(size_t row, size_t col, T init)
            : v(row, col, init), g(row, col, 0) {}

    template<tparams>
    DiffDenseMatrix::DenseMatrix(initializer_list list) : v(std::move(list)) {
        g = GradMatrix(v.getRow(), v.getCol(), 0);
    }

    template<tparams>
    template<Matrix M>
    DiffDenseMatrix::DenseMatrix(const M& mat) : DenseMatrix(mat.getRow(), mat.getCol()) {
        mat.assignTo(*this);
    }

    template<tparams>
    template<Vector V>
    DiffDenseMatrix::DenseMatrix(const V& vec) : DenseMatrix(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.assignTo(col);
    }

    template<tparams>
    template<RandomGenerator R>
    inline void DiffDenseMatrix::random_uniform() {
        *this = random_uniform<R>(getRow(), getCol());
    }

    template<tparams>
    template<RandomGenerator R>
    inline void DiffDenseMatrix::random_normal() {
        *this = random_normal<R>(getRow(), getCol());
    }

    template<tparams>
    template<class Distribution, RandomGenerator R>
    inline void DiffDenseMatrix::random_any(Distribution& dist) {
        *this = random_any<Distribution, R>(getRow(), getCol(), dist);
    }

    template<tparams>
    void DiffDenseMatrix::resize(size_t row, size_t col) {
        v.resize(row, col);
        g.resize(row, col);
    }

    template<tparams>
    void DiffDenseMatrix::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<tparams>
    void DiffDenseMatrix::swap_row(size_t r1, size_t r2) noexcept {
        v.swap_row(r1, r2);
        g.swap_row(r1, r2);
    }

    template<tparams>
    void DiffDenseMatrix::swap_col(size_t c1, size_t c2) noexcept {
        v.swap_col(c1, c2);
        g.swap_col(c1, c2);
    }

    template<tparams>
    inline DiffDenseMatrix::PtrTy DiffDenseMatrix::data_ptr(size_t row, size_t col) noexcept {
        assert(row < getRow() && "[Error]: Index out of range");
        assert(col < getCol() && "[Error]: Index out of range");
        return PtrTy(v.data_ptr(row, col), g.data_ptr(row, col));
    }

    template<tparams>
    inline DiffDenseMatrix::ConstPtrTy DiffDenseMatrix::data_ptr(size_t row, size_t col) const noexcept {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<tparams>
    DiffDenseMatrix DiffDenseMatrix::unitMatrix(size_t order) {
        This result(order, order);
        result.toUnitMatrix();
        return result;
    }

    template<tparams>
    template<RandomGenerator R>
    inline auto DiffDenseMatrix::random_uniform(size_t row, size_t col) {
        return This(ValueMatrix::template random_uniform<R>(row, col));
    }

    template<tparams>
    template<RandomGenerator R>
    inline auto DiffDenseMatrix::random_normal(size_t row, size_t col) {
        return This(ValueMatrix::template random_normal<R>(row, col));
    }

    template<tparams>
    template<class Distribution, RandomGenerator R>
    inline auto DiffDenseMatrix::random_any(size_t row, size_t col, Distribution& dist) {
        return This(ValueMatrix::template random_any<Distribution, R>(row, col, dist));
    }

#undef DiffDenseMatrix
#undef tparams
}
