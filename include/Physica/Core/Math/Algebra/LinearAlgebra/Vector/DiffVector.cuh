/*
 * Copyright 2024-2026 Weibo He.
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

#include "DenseVector.cuh"
#include "DiffVector.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order>
    class device_obj<DenseVector<Diff<T, Mode, Order>>>
            : public device_obj<CompactVector<DenseVector<Diff<T, Mode, Order>>>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<device_obj<DenseVector<Diff<T, Mode, Order>>>>, Empty>::type {
        using host_obj = DenseVector<Diff<T, Mode, Order>>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        template<Scalar U>
        using rebind_scalar = device_obj<DenseVector<U>>;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    private:
        using ValueVector = device_obj<VectorND<T>>;
        using GradType = ScalarType::GradType;
        using GradVector = std::conditional<Order == 1, ValueVector, device_obj<VectorND<GradType>>>::type;
        using initializer_list = std::initializer_list<T>;

        ValueVector v;
        GradVector g;
    public:
        device_obj() = default;
        explicit device_obj(size_t length);
        device_obj(size_t length, T init);
        device_obj(size_t length, ScalarType init) requires(isForwardDiff());
        template<Vector V>
        explicit(isReverseDiff()) device_obj(const V& v_) requires(!ReverseDiff<V>);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        void resize(this auto&, const Vector auto& x);
        void resize(this auto&, size_t size);

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;
        using Base::toHost;
        using Base::toHostAsync;

        template<RNG R>
        void random_uniform(this auto&);
        template<RNG R>
        void random_normal(this auto&);
        template<RNG R>
        void random_any(this auto&, auto& distribution);

        void zero_grad(this auto&);
        void junk(this auto&);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept;

        [[nodiscard]] __host__ __device__ auto&& values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto&& grads(this auto&&) noexcept;
        /* Static members */
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t len, auto& distribution);
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order>
    class Traits<device_obj<DenseVector<Diff<T, Mode, Order>>>> : public Traits<DenseVector<Diff<T, Mode, Order>>> {
        static_assert(!Diffable<T>, "[Error]: Nested Diff<> is not allowed");
    };
}

#include "DiffVectorImpl/DiffVectorImpl.cuh"
