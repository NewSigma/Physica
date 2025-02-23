/*
 * Copyright 2023-2025 Weibo He.
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

#include "RValueTensor.h"
#include "LValueTensorImpl/LTensorBlock.h"

namespace Physica {
    /**
     * Notes:
     * Right is positive direction of x.
     * Back is positive direction of y.
     * Top is positive direction of z.
     */
    template<class Derived>
    class LValueTensor : public RValueTensor<Derived> {
        using Base = RValueTensor<Derived>;
    public:
        using typename Base::ScalarType;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using RefTy = ScalarType::RefTy;
        using ConstRefTy = ScalarType::ConstRefTy;
    public:
        template<Tensor T>
        Derived& operator=(const T& other);
        Derived& operator=(const ScalarType& s);
        [[nodiscard]] ScalarType& operator()(size_t x, size_t y, size_t z) { return *data_ptr({x, y, z}); }
        [[nodiscard]] const ScalarType& operator()(size_t x, size_t y, size_t z) const { return *data_ptr({x, y, z}); }
        [[nodiscard]] ScalarType& operator()(Index3D index) { return *data_ptr(index); }
        [[nodiscard]] const ScalarType& operator()(Index3D index) const { return *data_ptr(index); }
        template<Scalar T> void operator+=(const T& s) { Base::getDerived() = Base::getDerived() + s; }
        template<Scalar T> void operator-=(const T& s) { Base::getDerived() = Base::getDerived() - s; }
        template<Scalar T> void operator*=(const T& s) { Base::getDerived() = Base::getDerived() * s; }
        template<Scalar T> void operator/=(const T& s) { Base::getDerived() = Base::getDerived() / s; }
        /* Operations */
        [[nodiscard]] inline LTensorBlock<Derived> block(Index3D from, Index3D count);
        [[nodiscard]] inline const LTensorBlock<Derived> block(Index3D from, Index3D count) const;

        [[nodiscard]] auto flatten() const { return FlattenTensor<Derived>(Base::getDerived()); }
        void resize(Index3D size) { Base::getDerived().resize(size); }
        template<RNG R> void random_uniform();
        template<RNG R> void random_normal();
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return operator()(index); }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(Index3D index) { return Base::getDerived().data_ptr(index); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(Index3D index) const { return Base::getDerived().data_ptr(index); }
        using Base::getDimX;
        using Base::getDimY;
        using Base::getDimZ;
        using Base::getDim;
        /* Static members */
        using Base::forPointInTensor;
        using Base::forPointIndexInTensor;
        using Base::forIndexInTensor;
    };

    template<Tensor T, Tensor U>
    inline void operator+=(LValueTensor<T>& g1, const U& g2) {
        g1.getDerived() = g1.getDerived() + g2.getDerived();
    }

    template<Tensor T, Tensor U>
    inline void operator-=(LValueTensor<T>& g1, const U& g2) {
        g1.getDerived() = g1.getDerived() - g2.getDerived();
    }
}

#include "LValueTensorImpl/LValueTensorImpl.h"
#include "LValueTensorImpl/FlattenTensor.h"
