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
    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::DenseMatrix(size_t row, size_t col) : values(row, col), grads(row, col) {}

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::DenseMatrix(size_t row, size_t col, T init)
            : values(row, col, init), grads(row, col, init) {}

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    template<Matrix M>
    DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::DenseMatrix(const M& mat) : DenseMatrix(mat.getRow(), mat.getCol()) {
        mat.assignTo(*this);
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    template<Vector V>
    DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::DenseMatrix(const V& vec) : DenseMatrix(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.assignTo(col);
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    void DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::resize(size_t row, size_t col) {
        values.resize(row, col);
        grads.resize(row, col);
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    void DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    void DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::swap_row(size_t r1, size_t r2) noexcept {
        values.swap_row(r1, r2);
        grads.swap_row(r1, r2);
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    void DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::swap_col(size_t c1, size_t c2) noexcept {
        values.swap_col(c1, c2);
        grads.swap_col(c1, c2);
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    inline DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::PtrTy
    DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::data_ptr(size_t row, size_t col) noexcept {
        return PtrTy(values.data_ptr(row, col), grads.data_ptr(row, col));
    }

    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    inline DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::ConstPtrTy
    DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col>::data_ptr(size_t row, size_t col) const noexcept {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
    ////////////////////////////////////////////////////////////////////////
    template<Scalar T, int Option, int Order>
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::Diff()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<Scalar T, int Option, int Order>
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::Diff(size_t row, size_t col)
            : traceSeg(TracerType::getInstance().pushSegment(row * col)) {}

    template<Scalar T, int Option, int Order>
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::Diff(PlainMatrix values)
            : traceSeg(TracerType::getInstance().pushSegment(values.flatten())) {}

    template<Scalar T, int Option, int Order>
    inline Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::ScalarType
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::calc(size_t row, size_t col) const {
        if constexpr (MatrixOption::isRowMatrix<This>())
            return traceSeg[row * getCol() + col];
        else
            return traceSeg[col * getRow() + row];
    }

    template<Scalar T, int Option, int Order>
    template<RandomGenerator R>
    inline void Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::random_uniform() {
        *this = random_uniform<R>(getRow(), getCol());
    }

    template<Scalar T, int Option, int Order>
    template<RandomGenerator R>
    inline void Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::random_normal() {
        *this = random_normal<R>(getRow(), getCol());
    }

    template<Scalar T, int Option, int Order>
    template<class Distribution, RandomGenerator R>
    inline void Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::random_any(Distribution& dist) {
        *this = random_any(getRow(), getCol(), dist);
    }

    template<Scalar T, int Option, int Order>
    template<RandomGenerator R>
    inline Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::random_uniform(
            size_t row, size_t col) {
        return This(PlainMatrix::template random_uniform<R>(row, col));
    }

    template<Scalar T, int Option, int Order>
    template<RandomGenerator R>
    inline Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::random_normal(
            size_t row, size_t col) {
        return This(PlainMatrix::template random_normal<R>(row, col));
    }

    template<Scalar T, int Option, int Order>
    template<class Distribution, RandomGenerator R>
    inline Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>
    Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>::random_any(
            size_t row, size_t col, Distribution& dist) {
        return This(PlainMatrix::random_any(row, col, dist));
    }
}
