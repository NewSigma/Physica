/*
 * Copyright 2026 Weibo He.
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

#include "../StridedVector.cuh"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<Vector V, size_t Length>
    class device_obj<StridedVectorBlock<V, Length>> : public device_obj<StridedVector<StridedVectorBlock<V, Length>>> {
        using host_obj = StridedVectorBlock<V, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedVector<host_obj>>;
        using Ref = add_device_obj<V>::type;
    private:
        Physica::PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(Ref vec_, size_t from, size_t to);
        __host__ __device__ device_obj(Ref vec_, size_t from);
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
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&&) noexcept;
    };

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<StridedVectorBlock<V, Length>>::device_obj(Ref vec_, size_t from, size_t to)
            : vec(asStruct(vec_)), from(from), to(to) {
        assert(from < to);
        assert(to <= vec.getDerived().getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<StridedVectorBlock<V, Length>>::device_obj(Ref vec_, size_t from)
            : device_obj(vec_, from, vec_.getLength()) {}

    template<Vector V, size_t Length>
    auto device_obj<StridedVectorBlock<V, Length>>::operator=(const This& obj) -> This& {
        Base::template operator=<This>(obj);
        return *this;
    }

    template<Vector V, size_t Length>
    auto device_obj<StridedVectorBlock<V, Length>>::operator=(This&& obj) noexcept -> This& {
        Base::template operator=<This>(std::move(obj));
        return *this;
    }

    template<Vector V, size_t Length>
    __host__ __device__ constexpr size_t device_obj<StridedVectorBlock<V, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector V, size_t Length>
    __host__ __device__ constexpr size_t device_obj<StridedVectorBlock<V, Length>>::getStride() const noexcept {
        return vec.getDerived().getStride();
    }

    template<Vector V, size_t Length>
    __host__ __device__ auto device_obj<StridedVectorBlock<V, Length>>::data_handle(this auto&& self) noexcept {
        return self.vec.getDerived().data_handle() + self.from * self.vec.getDerived().getStride();
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<device_obj<StridedVectorBlock<V, Length>>> : public Traits<StridedVectorBlock<V, Length>> {};
}
