/*
 * Copyright 2021-2026 Weibo He.
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

#include "TensorImpl/CompactTensor.h"

namespace Physica {
    template<Scalar T, int... Dims>
    class DenseTensor : public CompactTensor<DenseTensor<T, Dims...>> {
        using This = DenseTensor<T, Dims...>;
        using Base = CompactTensor<This>;
    public:
        using Base::NDim;
        using typename Base::IndexType;
    private:
        ArrayND<T, Dims...> storage;
    public:
        DenseTensor() = default;
        explicit DenseTensor(ArrayND<T, Dims...> storage) noexcept;
        DenseTensor(IndexType shape, auto&&... args);
        DenseTensor(std::integral auto... dims);
        DenseTensor(const Tensor auto& x);
        DenseTensor(const This&) = default;
        DenseTensor(This&&) noexcept = default;
        ~DenseTensor() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using Base::resize;
        void resize(this auto&, IndexType shape);
        void reserve(this auto&, size_t size) noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        using Base::random_normal;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto data_handle(this auto&&) noexcept;
        [[nodiscard]] IndexType getShape() const noexcept;
        [[nodiscard]] constexpr size_t dim(int index) const noexcept;
        [[nodiscard]] auto&& asArray(this auto&&) noexcept;
        [[nodiscard]] size_t getSize() const noexcept;
        /* Static members */
        template<RNG R>
        static DenseTensor random_uniform(IndexType shape);
        template<RNG R>
        static DenseTensor random_normal(IndexType shape);
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T, int... Dims>
    void swap(DenseTensor<T, Dims...>& x, DenseTensor<T, Dims...>& y) noexcept {
        x.swap(y);
    }

    template<Scalar T> using Tensor3D = DenseTensor<T, 3>;
    template<Scalar T> using Tensor4D = DenseTensor<T, 4>;
    template<Scalar T> using TensorND = DenseTensor<T>;
}

namespace Physica {
    template<Scalar T, int... Dims>
    class Traits<DenseTensor<T, Dims...>> {
        static_assert(sizeof...(Dims) > 0, "[Error]: Dims is not specified");
    public:
        using ScalarType = T;
        constexpr static int NDim = ArrayND<T, Dims...>::NDim;
    };
}

#include "TensorImpl/DenseTensorImpl.h"
