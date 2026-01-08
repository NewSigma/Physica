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

#include "../ContinuousVector.cuh"

namespace Physica {
    template<Vector T, size_t Length>
    class device_obj<ContinuousVectorBlock<T, Length>> : public device_obj<ContinuousVector<ContinuousVectorBlock<T, Length>>> {
        using host_obj = ContinuousVectorBlock<T, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;
    private:
        Physica::PlainStruct<device_obj<T>> vec;
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
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
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
    auto device_obj<ContinuousVectorBlock<T, Length>>::values(this auto&& self) noexcept {
        auto&& v = self.vec.getDerived().values();
        using V1 = remove_device_obj<std::remove_reference_t<decltype(v)>>::type;
        return device_obj<ContinuousVectorBlock<V1, Length>>(v, self.from, self.to);
    }

    template<Vector T, size_t Length>
    template<int GradOrder>
    auto device_obj<ContinuousVectorBlock<T, Length>>::grads(this auto&& self) noexcept {
        auto&& g = self.vec.template grads<GradOrder>();
        using V1 = remove_device_obj<std::remove_reference_t<decltype(g)>>::type;
        return ContinuousVectorBlock<V1, Length>(g, self.from, self.to);
    }

    template<Vector T, size_t Length>
    __host__ __device__ size_t device_obj<ContinuousVectorBlock<T, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }

    template<Vector T, size_t Length>
    __host__ __device__ auto device_obj<ContinuousVectorBlock<T, Length>>::data(this auto&& self) noexcept {
        return self.vec.getDerived().data() + self.from;
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<device_obj<ContinuousVectorBlock<T, Length>>> : public Traits<ContinuousVectorBlock<T, Length>> {};
}
