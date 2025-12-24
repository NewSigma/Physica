/*
 * Copyright 2020-2025 Weibo He.
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

#include "Physica/Core/Utils/Container/Array.h"
#include "VectorImpl/ContinuousVector.h"

namespace Physica {
    template<Scalar T, size_t Length = Dynamic, class Allocator = HostAllocator<T, alignof(typename BestPacket<T, Length>::Type)>>
    class DenseVector : public ContinuousVector<DenseVector<T, Length, Allocator>>
                      , public CRCoro<DenseVector<T, Length, Allocator>>
                      , public Array<T, Length, Allocator> {
        constexpr static size_t DefaultAlign = alignof(typename BestPacket<T, Length>::Type);
        static_assert(std::allocator_traits<Allocator>::Align % DefaultAlign == 0, "[Error]: Bad alignment for SIMD");
        using This = DenseVector<T, Length, Allocator>;
        using Coro = CRCoro<This>;
    public:
        using promise_type = Coro::promise_type;
        using Base = ContinuousVector<This>;
        using Storage = Array<T, Length, Allocator>;
        using device_obj_type = device_obj<This>;
        using typename Base::ScalarType;
        using Base::SizeAtCompile;
        using Base::isReverseDiff;
    protected:
        using typename Base::Trv;
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    public:
        DenseVector() = default;
        explicit DenseVector(size_t length);
        DenseVector(size_t length, const T& init);
        DenseVector(std::initializer_list<T> list);
        explicit DenseVector(Storage array) noexcept;
        DenseVector(const Vector auto& v);
        DenseVector(const This&) = default;
        DenseVector(This&&) noexcept = default;
        ~DenseVector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        void resize(const Vector auto& x);
        using Storage::resize;
        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        using Base::toDevice;
        using Base::toDeviceAsync;

        using Base::begin;
        using Base::cbegin;
        using Base::end;
        using Base::cend;
        using Base::rbegin;
        using Base::crbegin;
        using Base::rend;
        using Base::crend;

        [[nodiscard]] T skew() const;
        [[nodiscard]] T excess_kurt() const;
        [[nodiscard]] T kurt() const;

        using Base::zeros;
        using Base::random_uniform;
        using Base::random_normal;
        using Base::random_any;

        using Base::read;
        using Base::write;
        using Storage::swap;

        using Coro::get_return_object;
        using Coro::initial_suspend;
        using Coro::final_suspend;
        using Coro::return_value;
        using Coro::unhandled_exception;
        /* Getters */
        using Storage::getLength;
        using Storage::data;
        using Storage::data_ptr;
        /* Static members */
        [[nodiscard]] static This zeros(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_uniform(const This& v1, const This& v2);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t len, auto& distribution);
        [[nodiscard]] static This linspace(T from, T to, size_t count);
        [[nodiscard]] static This read_hdf5(const H5Loc& loc, const char* name);
        [[nodiscard]] static This read(size_t length, const T* __restrict p);
        [[nodiscard]] static This generate(std::invocable<size_t> auto fn);
        [[nodiscard]] static This generate(size_t length, std::invocable<size_t> auto fn);
    };

    template<Scalar T> using Vector1D = DenseVector<T, 1>;
    template<Scalar T> using Vector2D = DenseVector<T, 2>;
    template<Scalar T> using Vector3D = DenseVector<T, 3>;
    template<Scalar T> using Vector4D = DenseVector<T, 4>;
    template<Scalar T> using VectorND = DenseVector<T>;
}

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    class Traits<DenseVector<T, Length, Allocator>> {
        static_assert(!Diffable<T>, "[Error]: Use diffable vector instead");
    public:
        using ScalarType = T;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;

        using ElemType = T;
    };
}

namespace std {
    template<Physica::Scalar T, size_t Length>
    struct tuple_size<Physica::DenseVector<T, Length>> : public integral_constant<std::size_t, Length> {};

    template<size_t I, Physica::Scalar T, size_t Length>
    struct tuple_element<I, Physica::DenseVector<T, Length>> {
        using type = T;
    };

    template<Physica::Scalar T, size_t Length>
    void swap(Physica::DenseVector<T, Length>& __restrict v1, Physica::DenseVector<T, Length>& __restrict v2) noexcept {
        v1.swap(v2);
    }
}

#include "VectorImpl/DenseVectorImpl.h"
