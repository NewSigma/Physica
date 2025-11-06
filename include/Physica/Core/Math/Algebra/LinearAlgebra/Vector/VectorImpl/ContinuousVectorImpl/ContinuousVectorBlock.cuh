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

#include "../ContinuousVector.cuh"

namespace Physica {
    template<Vector T, size_t Length>
    class device_obj<ContinuousVectorBlock<T, Length>> : public device_obj<ContinuousVector<ContinuousVectorBlock<T, Length>>> {
        using host_obj = ContinuousVectorBlock<T, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;
        using DeviceVector = device_obj<T>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
        using typename Base::RefTy;
        using typename Base::ConstRefTy;
    private:
        Physica::PlainStruct<DeviceVector> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousVector<T>>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(device_obj<ContinuousVector<T>>& vec_, size_t from_);
        device_obj(const This& block) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& obj);
        This& operator=(This&& obj) noexcept;
        [[nodiscard]] __device__ RefTy operator[](size_t index);
        [[nodiscard]] __device__ ConstRefTy operator[](size_t index) const;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) noexcept { return vec.getDerived().data() + from + index; }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const noexcept { return vec.getDerived().data() + from + index; }
    };

    template<Vector T, size_t Length>
    __host__ __device__ device_obj<ContinuousVectorBlock<T, Length>>::device_obj(
            device_obj<ContinuousVector<T>>& vec_,
            size_t from_,
            size_t to_) : vec(asStruct(vec_.getDerived())), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getDerived().getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    __host__ __device__ device_obj<ContinuousVectorBlock<T, Length>>::device_obj(
            device_obj<ContinuousVector<T>>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}

    template<Vector T, size_t Length>
    auto device_obj<ContinuousVectorBlock<T, Length>>::operator=(const This& obj) -> This& {
        Base::template operator=<This>(obj);
        return *this;
    }

    template<Vector T, size_t Length>
    auto device_obj<ContinuousVectorBlock<T, Length>>::operator=(This&& obj) noexcept -> This& {
        Base::template operator=<This>(std::move(obj));
        return *this;
    }

    template<Vector T, size_t Length>
    __device__ auto device_obj<ContinuousVectorBlock<T, Length>>::operator[](size_t index) -> RefTy {
        assert((index + from) < to);
        return vec.getDerived()[index + from];
    }

    template<Vector T, size_t Length>
    __device__ auto device_obj<ContinuousVectorBlock<T, Length>>::operator[](size_t index) const -> ConstRefTy {
        assert((index + from) < to);
        return vec.getDerived()[index + from];
    }

    template<Vector T, size_t Length>
    __host__ __device__ size_t device_obj<ContinuousVectorBlock<T, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<device_obj<ContinuousVectorBlock<T, Length>>> : public Traits<ContinuousVectorBlock<T, Length>> {};
}
