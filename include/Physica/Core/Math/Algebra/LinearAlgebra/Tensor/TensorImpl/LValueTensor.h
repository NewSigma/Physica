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
    template<class Derived>
    class LValueTensor : public RValueTensor<Derived> {
        using This = LValueTensor<Derived>;
        using Base = RValueTensor<Derived>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using RefTy = ScalarType::RefTy;
        using ConstRefTy = ScalarType::ConstRefTy;
    public:
        ~LValueTensor() = default;
        /* Operators */
        template<Tensor T>
        Derived& operator=(const T& other);

        template<Scalar U>
        Derived& operator=(const U& x) requires(!isReverseDiff || !ReverseDiff<U>);
        template<Scalar U> void operator+=(const U& x) { Base::getDerived() = Base::getDerived() + x; }
        template<Scalar U> void operator-=(const U& x) { Base::getDerived() = Base::getDerived() - x; }
        template<Scalar U> void operator*=(const U& x) { Base::getDerived() = Base::getDerived() * x; }
        template<Scalar U> void operator/=(const U& x) { Base::getDerived() = Base::getDerived() / x; }

        template<Tensor X> inline void operator+=(const X& x);
        template<Tensor X> inline void operator-=(const X& x);

        [[nodiscard]] ScalarType& operator()(size_t x, size_t y, size_t z) { return *data_ptr({x, y, z}); }
        [[nodiscard]] const ScalarType& operator()(size_t x, size_t y, size_t z) const { return *data_ptr({x, y, z}); }
        [[nodiscard]] ScalarType& operator()(Index3D index) { return *data_ptr(index); }
        [[nodiscard]] const ScalarType& operator()(Index3D index) const { return *data_ptr(index); }
        /* Operations */
        [[nodiscard]] ScalarType calc(Index3D index) const { return operator()(index); }

        [[nodiscard]] inline LTensorBlock<Derived> block(Index3D from, Index3D count);
        [[nodiscard]] inline const LTensorBlock<Derived> block(Index3D from, Index3D count) const;

        [[nodiscard]] auto flatten();
        [[nodiscard]] const auto flatten() const;

        template<RNG R> void random_uniform();
        template<RNG R> void random_normal();
        /* Getters */
        [[nodiscard]] ScalarType* data_ptr(Index3D index) { return Base::getDerived().data_ptr(index); }
        [[nodiscard]] const ScalarType* data_ptr(Index3D index) const { return Base::getDerived().data_ptr(index); }
        /* Static members */
        using Base::forPointInTensor;
        using Base::forPointIndexInTensor;
        using Base::forIndexInTensor;
    protected:
        LValueTensor() = default;
        LValueTensor(const This&) = default;
        LValueTensor(This&&) noexcept = default;
    };
}

#include "LValueTensorImpl/LValueTensorImpl.h"
#include "LValueTensorImpl/Flatten.h"
