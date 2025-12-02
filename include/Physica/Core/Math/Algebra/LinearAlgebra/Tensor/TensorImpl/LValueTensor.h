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
        using typename Base::IndexType;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using PtrTy = T::PtrTy;
        using ConstPtrTy = T::ConstPtrTy;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
    public:
        ~LValueTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        template<Scalar U>
        Derived& operator=(const U& x) requires(!isReverseDiff || !ReverseDiff<U>);
        void operator+=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() + x; }
        void operator-=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() - x; }
        void operator*=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() * x; }
        void operator/=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() / x; }

        Derived& operator=(const Tensor auto& other);
        void operator+=(const Tensor auto& x);
        void operator-=(const Tensor auto& x);

        [[nodiscard]] RefTy& operator()(size_t x, size_t y, size_t z) { return *data_ptr({x, y, z}); }
        [[nodiscard]] ConstRefTy& operator()(size_t x, size_t y, size_t z) const { return *data_ptr({x, y, z}); }
        [[nodiscard]] RefTy& operator()(Index3D index) { return *data_ptr(index); }
        [[nodiscard]] ConstRefTy& operator()(Index3D index) const { return *data_ptr(index); }
        /* Operations */
        [[nodiscard]] T calc(Index3D index) const { return operator()(index); }

        void forND(std::invocable<T&, IndexType> auto fn);
        void forND(std::invocable<const T&, IndexType> auto fn) const;

        [[nodiscard]] auto fiber(int dim, IndexType index) noexcept;
        [[nodiscard]] const auto fiber(int dim, IndexType index) const noexcept;
        [[nodiscard]] auto slice(int dimRow, int dimCol, IndexType index) noexcept;
        [[nodiscard]] const auto slice(int dimRow, int dimCol, IndexType index) const noexcept;
        [[nodiscard]] LTensorBlock<Derived> block(Index3D from, Index3D count) noexcept;
        [[nodiscard]] const LTensorBlock<Derived> block(Index3D from, Index3D count) const noexcept;

        [[nodiscard]] auto flatten();
        [[nodiscard]] const auto flatten() const;

        template<RNG R> void random_uniform();
        template<RNG R> void random_normal();
        /* Getters */
        [[nodiscard]] PtrTy data_ptr(Index3D index) noexcept { return Base::getDerived().data_ptr(index); }
        [[nodiscard]] ConstPtrTy data_ptr(Index3D index) const noexcept { return Base::getDerived().data_ptr(index); }
        /* Static members */
        using Base::forPointIndexInTensor;
    protected:
        LValueTensor() = default;
        LValueTensor(const This&) = default;
        LValueTensor(This&&) noexcept = default;
    };
}

#include "LValueTensorImpl/LValueTensorImpl.h"
#include "LValueTensorImpl/Flatten.h"
