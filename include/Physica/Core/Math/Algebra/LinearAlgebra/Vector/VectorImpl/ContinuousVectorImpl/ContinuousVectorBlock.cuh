/*
 * Copyright 2023-2024 Weibo He.
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
    template<Vector T, size_t Length>
    class device_obj<ContinuousVectorBlock<T, Length>> : public device_obj<ContinuousVector<ContinuousVectorBlock<T, Length>>> {
        using host_obj = ContinuousVectorBlock<T, Length>;
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
        device_obj(const device_obj& block) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        inline device_obj& operator=(const device_obj& obj);
        inline device_obj& operator=(device_obj&& obj) noexcept;
        [[nodiscard]] __device__ inline RefTy operator[](size_t index);
        [[nodiscard]] __device__ inline ConstRefTy operator[](size_t index) const;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) { return vec.getDerived().data() + from + index; }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const { return vec.getDerived().data() + from + index; }
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
    inline device_obj<ContinuousVectorBlock<T, Length>>&
    device_obj<ContinuousVectorBlock<T, Length>>::operator=(const device_obj<ContinuousVectorBlock<T, Length>>& obj) {
        Base::operator=(obj);
        return *this;
    }
    
    template<Vector T, size_t Length>
    inline device_obj<ContinuousVectorBlock<T, Length>>&
    device_obj<ContinuousVectorBlock<T, Length>>::operator=(device_obj<ContinuousVectorBlock<T, Length>>&& obj) noexcept {
        Base::operator=(std::move(obj));
        return *this;
    }

    template<Vector T, size_t Length>
    __device__ inline device_obj<ContinuousVectorBlock<T, Length>>::RefTy
    device_obj<ContinuousVectorBlock<T, Length>>::operator[](size_t index) {
        assert((index + from) < to);
        return vec.getDerived()[index + from];
    }

    template<Vector T, size_t Length>
    __device__ inline device_obj<ContinuousVectorBlock<T, Length>>::ConstRefTy
    device_obj<ContinuousVectorBlock<T, Length>>::operator[](size_t index) const {
        assert((index + from) < to);
        return vec.getDerived()[index + from];
    }

    template<Vector T, size_t Length>
    template<Side Owner>
    __host__ __device__ inline size_t device_obj<ContinuousVectorBlock<T, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<Core::device_obj<Core::ContinuousVectorBlock<T, Length>>> : public Traits<Core::ContinuousVectorBlock<T, Length>> {};
}
