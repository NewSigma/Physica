/*
 * Copyright 2023-2026 Weibo He.
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

#include "../LValueVector.cuh"

namespace Physica {
    template<class Derived>
    template<Scalar U>
    __host__ __device__ device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const U& x) {
        Base::static_assert_assign(x);
        if (IsHost()) {
            constexpr int WarpSize = Physica::CUDADevAttr::WarpSize;
            const int numBlock = (Base::getLength() + WarpSize - 1) / WarpSize;
            const int numThread = WarpSize;
            auto func = [target_ = asStruct(Base::getDerived()), x] __device__() mutable {
                const unsigned int delta = gridDim.x * blockDim.x;
                const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
                auto& target = target_.getDerived();
                const size_t length = target.getLength();
                for (unsigned int shift = 0; shift < length; shift += delta) {
                    const unsigned int index = id + shift;
                    if (index < length)
                        target[index] = x;
                }
            };
            CUDAExecutor::launch<WarpSize>(func, KernelConfig(numBlock, numThread));
            return Base::getDerived();
        }
        else if constexpr (IsDevice()) {
            if constexpr (!std::same_as<T, U>)
                return operator=(T(x));
            else {
                for (size_t i = 0; i < Base::getLength(); ++i)
                    (*this)[i] = x;
                return Base::getDerived();
            }
        }
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator+=(const Scalar auto& x) {
        auto& v = Base::getDerived();
        (v + x).assign(v);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator-=(const Scalar auto& x) {
        auto& v = Base::getDerived();
        (v - x).assign(v);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator*=(const Scalar auto& x) {
        auto& v = Base::getDerived();
        (v * x).assign(v);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator/=(const Scalar auto& x) {
        auto& v = Base::getDerived();
        (v / x).assign(v);
    }

    template<class Derived>
    __host__ __device__ device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const Vector auto& v) {
        using V = std::remove_cvref_t<decltype(v)>;
        if constexpr (std::is_same<Derived, V>::value)
            assert(this != &v && "[Error]: Self assign is likely a bug");
        Base::assert_assign(v);

        auto& x = Base::getDerived();
        if constexpr (IsHost())
            x.resize(v);
        v.assign(x);
        return x;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator+=(const Vector auto& v) {
        auto& v0 = Base::getDerived();
        (v0 + v).assign(v0);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator-=(const Vector auto& v) {
        assert(Base::getLength() == v.getLength());
        Base::getDerived() += -v;
    }

    template<class Derived>
    __device__ decltype(auto) device_obj<LValueVector<Derived>>::operator[](this auto&& self, size_t index) {
        return *self.data_ptr(index);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename T::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        if constexpr (Scalar<U>) {
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad);
        }
        else {
            static_assert(Vector<U>, "[Error]: Unexpected type");
            assert(Base::getLength() == grad.getLength());
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad.calc(i));
        }
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ const auto device_obj<LValueVector<Derived>>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ const auto device_obj<LValueVector<Derived>>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ const auto device_obj<LValueVector<Derived>>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueVector<Derived>>::random_uniform() {
        Derived::template random_uniform<R>(Base::getLength()).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueVector<Derived>>::random_normal() {
        Derived::template random_normal<R>(Base::getLength()).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueVector<Derived>>::random_any(auto& distribution) {
        Derived::template random_any<R>(Base::getLength(), distribution).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength());
        return self.getDerived().data_ptr(index);
    }
}