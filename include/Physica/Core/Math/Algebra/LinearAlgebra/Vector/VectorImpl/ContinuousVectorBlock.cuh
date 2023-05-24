/*
 * Copyright 2023 WeiBo He.
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
    template<class VectorType, size_t Length>
    class device_obj<ContinuousVectorBlock<VectorType, Length>> : public device_obj<ContinuousVector<ContinuousVectorBlock<VectorType, Length>>> {
        using DeviceVector = device_obj<VectorType>;
    public:
        using host_obj = ContinuousVectorBlock<VectorType, Length>;
        using Base = device_obj<ContinuousVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        Physica::PlainStruct<DeviceVector> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousVector<VectorType>>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(device_obj<ContinuousVector<VectorType>>& vec_, size_t from_);
        device_obj(const device_obj& block) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        inline device_obj& operator=(const device_obj& obj);
        inline device_obj& operator=(device_obj&& obj) noexcept;
        [[nodiscard]] __device__ inline ScalarType& operator[](size_t index);
        [[nodiscard]] __device__ inline const ScalarType& operator[](size_t index) const;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept;
        [[nodiscard]] __device__ inline ScalarType* data();
        [[nodiscard]] __host__ __device__ inline const ScalarType* data() const;
    };

    template<class VectorType, size_t Length>
    __host__ __device__ device_obj<ContinuousVectorBlock<VectorType, Length>>::device_obj(
            device_obj<ContinuousVector<VectorType>>& vec_,
            size_t from_,
            size_t to_) : vec(asStruct(vec_.getDerived())), from(from_), to(to_) {}
    
    template<class VectorType, size_t Length>
    __host__ __device__ device_obj<ContinuousVectorBlock<VectorType, Length>>::device_obj(
            device_obj<ContinuousVector<VectorType>>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}

    template<class VectorType, size_t Length>
    inline device_obj<ContinuousVectorBlock<VectorType, Length>>&
    device_obj<ContinuousVectorBlock<VectorType, Length>>::operator=(const device_obj<ContinuousVectorBlock<VectorType, Length>>& obj) {
        Base::operator=(obj);
        return *this;
    }
    
    template<class VectorType, size_t Length>
    inline device_obj<ContinuousVectorBlock<VectorType, Length>>&
    device_obj<ContinuousVectorBlock<VectorType, Length>>::operator=(device_obj<ContinuousVectorBlock<VectorType, Length>>&& obj) noexcept {
        Base::operator=(std::move(obj));
        return *this;
    }

    template<class VectorType, size_t Length>
    __device__ inline typename device_obj<ContinuousVectorBlock<VectorType, Length>>::ScalarType&
    device_obj<ContinuousVectorBlock<VectorType, Length>>::operator[](size_t index) {
        assert((index + from) < to);
        return vec.getDerived()[index + from];
    }

    template<class VectorType, size_t Length>
    __device__ inline const typename device_obj<ContinuousVectorBlock<VectorType, Length>>::ScalarType&
    device_obj<ContinuousVectorBlock<VectorType, Length>>::operator[](size_t index) const {
        assert((index + from) < to);
        return vec.getDerived()[index + from];
    }

    template<class VectorType, size_t Length>
    __host__ __device__ inline size_t device_obj<ContinuousVectorBlock<VectorType, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }

    template<class VectorType, size_t Length>
    __device__ inline typename device_obj<ContinuousVectorBlock<VectorType, Length>>::ScalarType*
    device_obj<ContinuousVectorBlock<VectorType, Length>>::data() {
        return vec.getDerived().data() + from;
    }

    template<class VectorType, size_t Length>
    __host__ __device__ inline const typename device_obj<ContinuousVectorBlock<VectorType, Length>>::ScalarType*
    device_obj<ContinuousVectorBlock<VectorType, Length>>::data() const {
        return vec.getDerived().data() + from;
    }
}
