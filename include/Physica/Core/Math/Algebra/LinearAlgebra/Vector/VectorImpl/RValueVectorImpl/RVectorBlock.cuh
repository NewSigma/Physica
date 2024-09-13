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
    template<class VectorType, size_t Length>
    class device_obj<RVectorBlock<VectorType, Length>> : public device_obj<RValueVector<RVectorBlock<VectorType, Length>>> {
        using host_obj = RVectorBlock<VectorType, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        const VectorType& vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType>>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType>>& vec_, size_t from_);
        device_obj(const This& block) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return vec.calc(index + from); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept;
    };

    template<class VectorType, size_t Length>
    __host__ __device__ device_obj<RVectorBlock<VectorType, Length>>::device_obj(
            const device_obj<RValueVector<VectorType>>& vec_, size_t from_, size_t to_) : vec(asStruct(vec_.getDerived())), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<class VectorType, size_t Length>
    __host__ __device__ device_obj<RVectorBlock<VectorType, Length>>::device_obj(
            const device_obj<RValueVector<VectorType>>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}

    template<class VectorType, size_t Length>
    __host__ __device__ inline size_t device_obj<RVectorBlock<VectorType, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }
}

namespace Physica {
    template<class VectorType, size_t Length>
    class Traits<Core::device_obj<Core::RVectorBlock<VectorType, Length>>> : public Traits<Core::RVectorBlock<VectorType, Length>> {};
}
