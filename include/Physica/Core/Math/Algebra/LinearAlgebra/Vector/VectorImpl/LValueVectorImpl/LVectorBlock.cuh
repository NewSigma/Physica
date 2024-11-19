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
    template<Vector T, size_t Length>
    class device_obj<LVectorBlock<T, Length>> : public device_obj<LValueVector<LVectorBlock<T, Length>>> {
        using host_obj = LVectorBlock<T, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
        using DeviceVector = device_obj<T>;
    public:
        using ScalarType = typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        Physica::PlainStruct<DeviceVector> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(device_obj<LValueVector<T>>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(device_obj<LValueVector<T>>& vec_, size_t from_);
        LVectorBlock(const LVectorBlock& block) = delete;
        LVectorBlock(LVectorBlock&&) noexcept = delete;
        ~LVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v) { Base::operator=(static_cast<const RValueVector<This>&>(v)); return *this; }
        This& operator=(This&& v) noexcept { Base::operator=(static_cast<const RValueVector<This>&>(v)); return *this; }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept { return to - from; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    };

    template<Vector T, size_t Length>
    __host__ __device__ device_obj<LVectorBlock<T, Length>>::device_obj(
            device_obj<LValueVector<T>>& vec_, size_t from_, size_t to_) : vec(asStruct(vec_)), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    __host__ __device__ device_obj<LVectorBlock<T, Length>>::device_obj(
            device_obj<LValueVector<T>>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}
    
    template<Vector T, size_t Length>
    __host__ __device__ inline size_t device_obj<LVectorBlock<T, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }

    template<Vector T, size_t Length>
    __host__ __device__ inline typename device_obj<LVectorBlock<T, Length>>::PtrTy
    device_obj<LVectorBlock<T, Length>>::data_ptr(size_t index) {
        assert((index + from) < to);
        return vec.getDerived().data_ptr(index);
    }

    template<Vector T, size_t Length>
    __host__ __device__ inline const typename device_obj<LVectorBlock<T, Length>>::ConstPtrTy
    device_obj<LVectorBlock<T, Length>>::data_ptr(size_t index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<Core::device_obj<Core::LVectorBlock<T, Length>>> : public Traits<Core::LVectorBlock<T, Length>> {};
}
