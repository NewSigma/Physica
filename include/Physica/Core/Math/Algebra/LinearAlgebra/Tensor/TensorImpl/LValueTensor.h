/*
 * Copyright 2023-2026 Weibo He.
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
    public:
        ~LValueTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        auto operator=(Scalar auto x) noexcept -> Derived&;
        void operator+=(Scalar auto x) noexcept;
        void operator-=(Scalar auto x) noexcept;
        void operator*=(Scalar auto x) noexcept;
        void operator/=(Scalar auto x) noexcept;

        auto operator=(const Tensor auto& other) -> Derived&;
        void operator+=(const Tensor auto& x);
        void operator-=(const Tensor auto& x);

        [[nodiscard]] decltype(auto) operator[](this auto&&, size_t x, size_t y, size_t z);
        [[nodiscard]] decltype(auto) operator[](this auto&&, Index3D index);
        /* Operations */
        [[nodiscard]] decltype(auto) calc(Index3D index) const { return operator[](index); }

        void forND(std::invocable<T&, IndexType> auto fn);
        void forND(std::invocable<const T&, IndexType> auto fn) const;

        [[nodiscard]] auto fiber(this auto&&, IndexVar auto...) noexcept;
        [[nodiscard]] auto slice(this auto&&, IndexVar auto...) noexcept;
        [[nodiscard]] auto block(this auto&&, Index3D from, Index3D count) noexcept;

        [[nodiscard]] auto flatten(this auto&&);

        void zero_grad() noexcept;
        template<RNG R> void random_uniform();
        template<RNG R> void random_normal();
        /* Getters */
        [[nodiscard]] auto data_ptr(this auto&&, Index3D index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isLValueTensor() noexcept { return true; }
    protected:
        LValueTensor() = default;
        LValueTensor(const This&) = default;
        LValueTensor(This&&) noexcept = default;
    };
}

#include "LValueTensorImpl/LValueTensorImpl.h"
#include "LValueTensorImpl/Flatten.h"
