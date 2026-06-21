/*
 * Copyright 2021-2026 Weibo He.
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

#include "../DenseMatrix.h"

namespace Physica {
    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(Storage storage) noexcept : storage(std::move(storage)) {}

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(size_t order) : DenseMatrix(order, order) {}

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(size_t row, size_t col)
            : storage(row, col) {}

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(size_t row, size_t col, T value)
            : storage(row, col, std::move(value)) {}

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(std::initializer_list<T> list) : storage(list) {}

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(std::initializer_list<VectorIniter> list)
            : DenseMatrix(MatrixMajor::isColMatrix<This>() ? list.begin()->getLength() : list.size(),
                          MatrixMajor::isColMatrix<This>() ? list.size() : list.begin()->getLength()) {
        size_t major = 0;
        for (auto& v : list) {
            assert(v.getLength() == Base::getMaxMinor());
            if constexpr (MatrixMajor::isColMatrix<This>())
                this->col(major) = v;
            else
                this->row(major) = v;
            major += 1;
        }
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(const Matrix auto& mat)
            : DenseMatrix(mat.getRow(), mat.getCol()) {
        mat.assign(*this);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    DenseMatrix<T, Major, Row, Col, Allocator>::DenseMatrix(const Vector auto& vec) : DenseMatrix(vec.getLength(), 1) {
        auto col = this->col(0);
        vec.assign(col);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::resize(size_t order) {
        storage.resize(order);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::resize(const Matrix auto& m, auto&&... args) {
        Base::resize(m, std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::resize(size_t row, size_t col, auto&&... args) {
        storage.resize(row, col, std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::reserve(size_t size) noexcept {
        storage.reserve(size);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto&& DenseMatrix<T, Major, Row, Col, Allocator>::flatten(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.storage).asArray();
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::zeros() noexcept {
        storage.zeros();
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::swap_row(size_t r1, size_t r2) noexcept {
        storage.swap_row(r1, r2);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::swap_col(size_t c1, size_t c2) noexcept {
        storage.swap_col(c1, c2);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void DenseMatrix<T, Major, Row, Col, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        storage.swap(obj.storage);
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::data(this auto&& self) noexcept {
        return self.storage.data();
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto&& DenseMatrix<T, Major, Row, Col, Allocator>::asArray(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.storage).asArray();
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    __host__ __device__ consteval size_t DenseMatrix<T, Major, Row, Col, Allocator>::getRowAtCompile() noexcept {
        return Row;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    __host__ __device__ consteval size_t DenseMatrix<T, Major, Row, Col, Allocator>::getColAtCompile() noexcept {
        return Col;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    __host__ __device__ consteval int DenseMatrix<T, Major, Row, Col, Allocator>::getMajor() noexcept {
        return Major;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::zeros(size_t row, size_t col) -> This {
        DenseMatrix result(row, col);
        result.zeros();
        return result;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::identity(size_t order) -> This {
        DenseMatrix result(order, order);
        result.toIdentity();
        return result;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    template<RNG R>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::random_uniform(size_t row, size_t col) -> This {
        This result(row, col);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    template<RNG R>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::random_normal(size_t row, size_t col) -> This {
        This result(row, col);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    template<RNG R>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::random_any(size_t row, size_t col, auto& distribution) -> This {
        This result(row, col);
        result.template random_any<R>(distribution);
        return result;
    }

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::meshgrid(const Vector auto& vecX, const Vector auto& vecY) -> std::pair<This, This> {
        using MatrixType = DenseMatrix<T, Major, Row, Col, Allocator>;
        const size_t row = vecY.getLength();
        const size_t col = vecX.getLength();
        MatrixType x(row, col);
        MatrixType y(row, col);
        for (size_t i = 0; i < x.getMaxMajor(); ++i) {
            for (size_t j = 0; j < x.getMaxMinor(); ++j) {
                x.refFromMajorMinor(i, j) = vecX.calc(MatrixMajor::colFromMajorMinor<MatrixType>(i, j));
                y.refFromMajorMinor(i, j) = vecY.calc(MatrixMajor::rowFromMajorMinor<MatrixType>(i, j));
            }
        }
        return std::make_pair(std::move(x), std::move(y));
    }
    /**
     * Helper function that communicates with C libraries.
     */
    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    auto DenseMatrix<T, Major, Row, Col, Allocator>::read(size_t row, size_t col, const T* __restrict p) noexcept -> This {
        return This(Storage::read(row, col, p));
    }
}
