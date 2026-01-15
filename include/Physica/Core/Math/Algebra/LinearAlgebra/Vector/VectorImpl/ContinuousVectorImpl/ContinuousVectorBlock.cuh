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
    template<Vector V, size_t Length>
    class device_obj<ContinuousVectorBlock<V, Length>> : public device_obj<ContinuousVector<ContinuousVectorBlock<V, Length>>> {
        using host_obj = ContinuousVectorBlock<V, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;
    private:
        Physica::PlainStruct<device_obj<V>> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(device_obj<V>& vec, size_t from, size_t to);
        __host__ __device__ device_obj(device_obj<V>& vec, size_t from);
        device_obj(const This& block) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& obj);
        This& operator=(This&& obj) noexcept;
        /* Operations */
        template<size_t Length_ = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(this auto&&, size_t to = Length_) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(this auto&&, size_t from, size_t to) noexcept;

        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
    };

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<ContinuousVectorBlock<V, Length>>::device_obj(device_obj<V>& vec, size_t from, size_t to)
            : vec(asStruct(vec)), from(from), to(to) {
        assert(from < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<ContinuousVectorBlock<V, Length>>::device_obj(device_obj<V>& vec, size_t from)
            : device_obj(vec, from, vec.getLength()) {}

    template<Vector V, size_t Length>
    auto device_obj<ContinuousVectorBlock<V, Length>>::operator=(const This& obj) -> This& {
        Base::template operator=<This>(obj);
        return *this;
    }

    template<Vector V, size_t Length>
    auto device_obj<ContinuousVectorBlock<V, Length>>::operator=(This&& obj) noexcept -> This& {
        Base::template operator=<This>(std::move(obj));
        return *this;
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    __host__ __device__ auto device_obj<ContinuousVectorBlock<V, Length>>::head(this auto&& self, size_t to) noexcept {
        return device_obj<ContinuousVectorBlock<V, Length_>>(self.vec, self.from, self.from + to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    __host__ __device__ auto device_obj<ContinuousVectorBlock<V, Length>>::tail(this auto&& self, size_t from) noexcept {
        return device_obj<ContinuousVectorBlock<V, Length_>>(self.vec, self.from + from, self.to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    __host__ __device__ auto device_obj<ContinuousVectorBlock<V, Length>>::segment(this auto&& self, size_t from, size_t to) noexcept {
        return device_obj<ContinuousVectorBlock<V, Length_>>(self.vec, self.from + from, self.from + to);
    }

    template<Vector V, size_t Length>
    auto device_obj<ContinuousVectorBlock<V, Length>>::values(this auto&& self) noexcept {
        auto&& v = self.vec.getDerived().values();
        using V1 = remove_device_obj<std::remove_reference_t<decltype(v)>>::type;
        return device_obj<ContinuousVectorBlock<V1, Length>>(v, self.from, self.to);
    }

    template<Vector V, size_t Length>
    template<int GradOrder>
    auto device_obj<ContinuousVectorBlock<V, Length>>::grads(this auto&& self) noexcept {
        auto&& g = self.vec.template grads<GradOrder>();
        using V1 = remove_device_obj<std::remove_reference_t<decltype(g)>>::type;
        return ContinuousVectorBlock<V1, Length>(g, self.from, self.to);
    }

    template<Vector V, size_t Length>
    __host__ __device__ size_t device_obj<ContinuousVectorBlock<V, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector V, size_t Length>
    __host__ __device__ auto device_obj<ContinuousVectorBlock<V, Length>>::data(this auto&& self) noexcept {
        return self.vec.getDerived().data() + self.from;
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<device_obj<ContinuousVectorBlock<V, Length>>> : public Traits<ContinuousVectorBlock<V, Length>> {};
}
