/*
 * Copyright 2023-2025 Weibo He.
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
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::operator=(const This& obj) -> This& {
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::operator=(This&& obj) -> This& {
        return *this = obj;
    }

    template<class Derived>
    __host__ __device__ device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const Vector auto& v) requires(CUDA<decltype(v)>) {
        using V = std::remove_cvref_t<decltype(v)>;
        if constexpr (std::is_same<Derived, V>::value)
            assert(this != &v && "[Error]: Self assign is likely a bug");
        auto& x = Base::getDerived();
        if constexpr (IsHost())
            x.resize(v.getLength());
        v.assign(x);
        return x;
    }

    template<class Derived>
    template<Scalar T>
    inline device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const T& x) {
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

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator+=(const Scalar auto& x) {
        Base::getDerived() = Base::getDerived() + x;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator-=(const Scalar auto& x) {
        Base::getDerived() = Base::getDerived() - x;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator*=(const Scalar auto& x) {
        Base::getDerived() = Base::getDerived() * x;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::operator/=(const Scalar auto& x) {
        Base::getDerived() = Base::getDerived() / x;
    }

    template<class Derived>
    __host__ __device__ inline void device_obj<LValueVector<Derived>>::operator+=(const Vector auto& v) requires(CUDA<decltype(v)>) {
        assert(Base::getLength() == v.getLength());
        Base::getDerived() = Base::getDerived() + v;
    }

    template<class Derived>
    __host__ __device__ inline void device_obj<LValueVector<Derived>>::operator-=(const Vector auto& v) requires(CUDA<decltype(v)>) {
        assert(Base::getLength() == v.getLength());
        Base::getDerived() += -v;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueVector<Derived>>::reverse(const auto& grad) const noexcept requires(isReverseDiff) {
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename ScalarType::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
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
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const auto device_obj<LValueVector<Derived>>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const auto device_obj<LValueVector<Derived>>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const auto device_obj<LValueVector<Derived>>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::data_ptr(size_t index) -> PtrTy {
        assert(index < Base::getLength());
        return Base::getDerived().data_ptr(index);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueVector<Derived>>::data_ptr(size_t index) const -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(index);
    }
}