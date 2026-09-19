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

#include "Physica/Core/Utils/Container/ArrayND.cuh"
#include "TensorImpl/LValueTensor.cuh"
#include "DenseTensor.h"

namespace Physica {
    template<Scalar T, int... Dims>
    class device_obj<DenseTensor<T, Dims...>>
            : public device_obj<LValueTensor<DenseTensor<T, Dims...>>>
            , public CRCoro<device_obj<DenseTensor<T, Dims...>>> {
        static_assert(!Diffable<T>, "[Error]: Use diffable tensor instead");
        using host_obj = DenseTensor<T, Dims...>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueTensor<host_obj>>;
        using Coro = CRCoro<This>;
        using Storage = device_obj<ArrayND<T, Dims...>>;
    public:
        using typename Base::IndexType;
        template<Scalar U>
        using rebind_scalar = device_obj<DenseTensor<U, Dims...>>;
    private:
        Storage storage;
    public:
        device_obj() = default;
        explicit device_obj(Storage storage) noexcept;
        __host__ __device__ device_obj(IndexType shape, auto&&... args);
        __host__ __device__ device_obj(std::integral auto... dims);
        device_obj(const Tensor auto& x);
        device_obj(const host_obj& obj);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(this auto&, IndexType shape);
        void reserve(this auto&, size_t size);

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        using Base::random_normal;
        using Base::random_uniform;

        void zeros(this auto&) noexcept;
        void junk(this auto&) noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t dim(int index) const noexcept { return storage.dim(index); }
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept { return storage.getShape(); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, const IndexType& indices) noexcept { return self.storage.data_ptr(indices); }
        [[nodiscard]] __host__ __device__ auto&& asArray(this auto&& self) noexcept { return self.storage.asArray(); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return storage.getSize(); }
        /* Static members */
        template<RNG R>
        [[nodiscard]] static This random_uniform(IndexType shape);
        template<RNG R>
        [[nodiscard]] static This random_normal(IndexType shape);
        /* Friends */
        friend class DenseTensor<T, Dims...>;
    };

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::toDeviceAsync() const {
        device_obj<This> result;
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::toDevice(device_obj<This>& obj) const {
        storage.toDevice(obj.storage);
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::toDeviceAsync(device_obj<This>& obj) const {
        storage.toDeviceAsync(obj.storage);
    }
}

namespace Physica {
    template<Scalar T, int... Dims>
    class Traits<device_obj<DenseTensor<T, Dims...>>> : public Traits<DenseTensor<T, Dims...>> {};
}

#include "TensorImpl/DenseTensorImpl.cuh"
