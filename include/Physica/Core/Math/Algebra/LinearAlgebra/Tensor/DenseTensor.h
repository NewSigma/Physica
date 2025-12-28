/*
 * Copyright 2021-2025 Weibo He.
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

#include "TensorImpl/LValueTensor.h"

namespace Physica {
    template<Scalar T, int... Dims>
    class DenseTensor : public LValueTensor<DenseTensor<T, Dims...>>, private ArrayND<T, Dims...> {
        using This = DenseTensor<T, Dims...>;
        using Base = LValueTensor<This>;
        using Storage = ArrayND<T, Dims...>;
    public:
        using Base::NDim;
        using typename Base::IndexType;
    public:
        DenseTensor() = default;
        DenseTensor(std::integral auto... dims);
        DenseTensor(IndexType shape, auto&&... args);
        DenseTensor(const This&) = default;
        DenseTensor(This&&) noexcept = default;
        ~DenseTensor() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Storage::operator[];
        /* Operations */
        using Storage::toIndex1D;
        using Storage::toIndexND;
        using Storage::forND;

        using Base::resize;
        void resize(IndexType shape);

        using Base::random_normal;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Storage::dim;
        using Storage::getShape;
        using Storage::data_ptr;
        using Storage::asArray;
        using Storage::getSize;
        /* Static members */
        template<RNG R>
        static DenseTensor random_uniform(IndexType shape);
        template<RNG R>
        static DenseTensor random_normal(IndexType shape);
    };

    template<Scalar T, int... Dims>
    void swap(DenseTensor<T, Dims...>& __restrict x, DenseTensor<T, Dims...>& __restrict y) noexcept {
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
