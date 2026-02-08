/*
 * Copyright 2025-2026 Weibo He.
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

#include "../DiffDenseMatrix.cuh"

namespace Physica {
#define tparams Scalar T, DiffMode Mode, int Order, int Major
#define DenseMatrix DenseMatrix<Diff<T, Mode, Order>, Major>

    template<tparams>
    device_obj<DenseMatrix>::device_obj(size_t row, size_t col) : v(row, col), g(row, col) {}

    template<tparams>
    device_obj<DenseMatrix>::device_obj(size_t row, size_t col, T init) : v(row, col, init), g(row, col, 0) {}

    template<tparams>
    device_obj<DenseMatrix>::device_obj(ValueMatrix v_, GradMatrix g_) : v(std::move(v_)), g(std::move(g_)) {}

    template<tparams>
    template<Vector V>
    device_obj<DenseMatrix>::device_obj(const V& vec) requires(!ReverseDiff<V>)  : device_obj(vec.getLength(), 1) {
        if constexpr (!Diffable<V>) {
            auto col = v.col(0);
            vec.assign(col);
            g.zeros();
        }
        else {
            auto col = this->col(0);
            vec.assign(col);
        }
    }

    template<tparams>
    template<Matrix M>
    device_obj<DenseMatrix>::device_obj(const M& mat) requires(!ReverseDiff<M>) : device_obj(mat.getRow(), mat.getCol()) {
        if constexpr (!Diffable<M>) {
            mat.assign(v);
            g.zeros();
        }
        else
            mat.assign(*this);
    }

    template<tparams>
    void device_obj<DenseMatrix>::zero_grad() {
        g.zeros();
    }

    template<tparams>
    __host__ __device__ void device_obj<DenseMatrix>::resize(size_t row, size_t col) {
        v.resize(row, col);
        g.resize(row, col);
    }

    template<tparams>
    void device_obj<DenseMatrix>::zeros() {
        v.zeros();
        g.zeros();
    }

    template<tparams>
    auto device_obj<DenseMatrix>::toHost() const -> host_obj {
        host_obj result = toHostAsync();
        CUDAContext::getInstance().wait();
        return result;
    }

    template<tparams>
    auto device_obj<DenseMatrix>::toHostAsync() const -> host_obj {
        host_obj result(getRow(), getCol());
        toHostAsync(result);
        return result;
    }

    template<tparams>
    void device_obj<DenseMatrix>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<tparams>
    void device_obj<DenseMatrix>::toHostAsync(host_obj& obj) const {
        v.toHostAsync(obj.v);
        g.toHostAsync(obj.g);
    }

    template<tparams>
    template<RNG R>
    void device_obj<DenseMatrix>::random_uniform() {
        *this = random_uniform<R>(getRow(), getCol());
    }

    template<tparams>
    template<RNG R>
    void device_obj<DenseMatrix>::random_normal() {
        *this = random_normal<R>(getRow(), getCol());
    }

    template<tparams>
    template<RNG R>
    void device_obj<DenseMatrix>::random_any(auto& distribution) {
        *this = random_any<R>(getRow(), getCol(), distribution);
    }

    template<tparams>
    void device_obj<DenseMatrix>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<tparams>
    __host__ __device__ auto device_obj<DenseMatrix>::data(this auto&& self) noexcept {
        constexpr bool IsConst = std::is_const<std::remove_reference_t<decltype(self)>>::value;
        using U = Diff<T, Mode, Order>;
        using RetTy = std::conditional<IsConst, typename U::ConstPtrTy, typename U::PtrTy>::type;
        return RetTy(self.v.data(), self.g.data());
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix>::random_uniform(size_t row, size_t col) -> This {
        return This(ValueMatrix::template random_uniform<R>(row, col));
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix>::random_normal(size_t row, size_t col) -> This {
        return This(ValueMatrix::template random_normal<R>(row, col));
    }

    template<tparams>
    template<RNG R>
    auto device_obj<DenseMatrix>::random_any(size_t row, size_t col, auto& distribution) -> This {
        return This(ValueMatrix::template random_any<R>(row, col, distribution));
    }

#undef tparams
#undef DenseMatrix

    template<Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col>
    auto DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>::toDevice() const {
        auto result = toDeviceAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col>
    auto DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>::toDeviceAsync() const {
        device_obj<This> result{};
        toDeviceAsync(result);
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col>
    void DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>::toDevice(device_obj<This>& obj) const {
        toDeviceAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col>
    void DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>::toDeviceAsync(device_obj<This>& obj) const {
        v.toDeviceAsync(obj.values());
        g.toDeviceAsync(obj.grads());
    }
}
