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
    template<Scalar T, int Dim = Dynamic>
    class DenseTensor : public LValueTensor<DenseTensor<T, Dim>>, private ArrayND<T, Dim> {
        using This = DenseTensor<T, Dim>;
        using Base = LValueTensor<This>;
        using Storage = ArrayND<T, Dim>;
    protected:
        using typename Base::IndexArray;
    public:
        DenseTensor() = default;
        DenseTensor(IndexArray shape, auto&&... args);
        DenseTensor(const This&) = default;
        DenseTensor(This&&) noexcept = default;
        ~DenseTensor() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Storage::operator();
        /* Operations */
        using Base::random_normal;
        void resize(IndexArray shape, auto&&... args);

        using Storage::toIndex1D;
        using Storage::toIndexND;
        using Storage::forND;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Storage::asArray;
        using Storage::data_ptr;
        using Storage::getShape;
        using Storage::getDim;
        using Storage::getSize;
        /* Static members */
        template<RNG R>
        static DenseTensor random_uniform(IndexArray shape);
        template<RNG R>
        static DenseTensor random_normal(IndexArray shape);
    };

    template<Scalar T> using Tensor3D = DenseTensor<T, 3>;
    template<Scalar T> using Tensor4D = DenseTensor<T, 4>;
    template<Scalar T> using TensorND = DenseTensor<T>;
}

namespace Physica {
    template<Scalar T, int Dim_>
    class Traits<DenseTensor<T, Dim_>> {
    public:
        using ScalarType = T;
        constexpr static int Dim = Dim_;
    };
}

namespace std {
    template<Physica::Scalar T, size_t Dim>
    void swap(Physica::DenseTensor<T, Dim>& __restrict x, Physica::DenseTensor<T, Dim>& __restrict y) noexcept {
        x.swap(y);
    }
}

#include "TensorImpl/DenseTensorImpl.h"
