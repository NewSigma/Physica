/*
 * Copyright 2024-2026 Weibo He.
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

namespace Physica {
#define tparams Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col
#define DiffDenseMatrix DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>

    template<tparams>
    DiffDenseMatrix::DenseMatrix(size_t order) : DenseMatrix(order, order) {}

    template<tparams>
    DiffDenseMatrix::DenseMatrix(size_t row, size_t col) : v(row, col) {
        if constexpr (isForwardDiff())
            g.resize(row, col);
    }

    template<tparams>
    DiffDenseMatrix::DenseMatrix(size_t row, size_t col, T init) : v(row, col, init), g(GradMatrix::zeros(row, col)) {}

    template<tparams>
    DiffDenseMatrix::DenseMatrix(size_t row, size_t col, ScalarType init) requires(isForwardDiff())
            : v(row, col, init.value()), g(row, col, init.grad()) {}

    template<tparams>
    DiffDenseMatrix::DenseMatrix(std::initializer_list<Tv> list) : v(list) {
        if constexpr (isForwardDiff())
            g.resize(v.getRow(), v.getCol());
        zero_grad();
    }

    template<tparams>
    DiffDenseMatrix::DenseMatrix(std::initializer_list<VectorIniter> list) : v(list) {
        if constexpr (isForwardDiff())
            g.resize(v.getRow(), v.getCol());
        zero_grad();
    }

    template<tparams>
    DiffDenseMatrix::DenseMatrix(ValueMatrix v_, GradMatrix g_) : v(std::move(v_)), g(std::move(g_)) {}

    template<tparams>
    DiffDenseMatrix::DenseMatrix(const Vector auto& vec) : DenseMatrix(vec.getLength(), 1) {
        using V = decltype(vec);
        static_assert(!ReverseDiff<V>, "[Error]: Copying a reverse diff object breaks the compute graph");
        if constexpr (!Diffable<V>) {
            vec.assign(v.col(0));
            zero_grad();
        }
        else
            vec.assign(this->col(0));
    }

    template<tparams>
    DiffDenseMatrix::DenseMatrix(const Matrix auto& mat) : DenseMatrix(mat.getRow(), mat.getCol()) {
        using M = decltype(mat);
        static_assert(!ReverseDiff<M>, "[Error]: Copying a reverse diff object breaks the compute graph");
        if constexpr (!Diffable<M>) {
            mat.assign(v);
            zero_grad();
        }
        else
            mat.assign(*this);
    }

    template<tparams>
    void DiffDenseMatrix::resize(size_t row, size_t col) {
        v.resize(row, col);
        if constexpr (isForwardDiff())
            g.resize(row, col);
        else
            assert(g.empty() && "[Error]: Reject resize-after-reverse");
    }

    template<tparams>
    void DiffDenseMatrix::reserve(size_t size) noexcept {
        v.reserve(size);
        if constexpr (isForwardDiff())
            g.reserve(size);
    }

    template<tparams>
    template<RNG R>
    void DiffDenseMatrix::random_uniform() {
        v.template random_uniform<R>();
        zero_grad();
    }

    template<tparams>
    template<RNG R>
    void DiffDenseMatrix::random_normal() {
        v.template random_normal<R>();
        zero_grad();
    }

    template<tparams>
    template<RNG R>
    void DiffDenseMatrix::random_any(auto& distribution) {
        v.template random_any<R>(distribution);
        zero_grad();
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
    void DiffDenseMatrix::zero_grad() {
        g.zeros();
    }

    template<tparams>
    auto DiffDenseMatrix::data(this auto&& self) noexcept {
        constexpr bool IsConst = std::is_const<std::remove_reference_t<decltype(self)>>::value;
        using U = Diff<T, Mode, Order>;
        using RetTy = std::conditional<IsConst, typename U::ConstPtrTy, typename U::PtrTy>::type;
        if constexpr (isReverseDiff())
            assert(!self.g.empty() && "[Error]: Grad is not ready");
        return RetTy(self.values().data(), self.grads().data());
    }

    template<tparams>
    auto&& DiffDenseMatrix::values(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.v);
    }

    template<tparams>
    template<int GradOrder>
    auto&& DiffDenseMatrix::grads(this auto&& self) noexcept {
        if constexpr (isReverseDiff()) {
            if (self.g.empty()) {
                auto& g1 = const_cast<GradMatrix&>(self.g);
                g1.resize(self.v);
                g1.zeros();
            }
        }

        if constexpr (GradOrder == 1)
            return forward_like<decltype(self)>(self.g);
        else
            return forward_like<decltype(self)>(self.g.template grads<GradOrder - 1>());
    }

    template<tparams>
    DiffDenseMatrix DiffDenseMatrix::identity(size_t order) {
        This result(order, order);
        result.toIdentity();
        return result;
    }

    template<tparams>
    template<RNG R>
    auto DiffDenseMatrix::random_uniform(size_t row, size_t col) {
        return This(ValueMatrix::template random_uniform<R>(row, col));
    }

    template<tparams>
    template<RNG R>
    auto DiffDenseMatrix::random_normal(size_t row, size_t col) {
        return This(ValueMatrix::template random_normal<R>(row, col));
    }

    template<tparams>
    template<RNG R>
    auto DiffDenseMatrix::random_any(size_t row, size_t col, auto& distribution) {
        return This(ValueMatrix::template random_any<R>(row, col, distribution));
    }

#undef DiffDenseMatrix
#undef tparams
}
