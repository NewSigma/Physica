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

#include "Physica/CRTPBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    template<class Derived> class LValueTensor;
    template<class> class FlattenL;
    template<class T> class RealTensor;
    template<class T> class ImagTensor;
    template<class T> class SquaredNormTensor;
    template<class T> class NormTensor;
    template<class T> class ValueTensor;    
    template<class T, int GradOrder> class GradTensor;

    template<class Derived>
    class RValueTensor : public CRTPBase<RValueTensor<Derived>> {
        static_assert(!DeviceObj<Derived>, "[Error]: device_obj<> must be outside RValueTensor<>");
        using This = RValueTensor<Derived>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
        constexpr static int NDim = Traits<Derived>::NDim;

        using IndexType = Array<size_t, NDim>;
    protected:
        using T = ScalarType;
    public:
        ~RValueTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Tensor auto& x) const;
        void assert_assign(const Tensor auto& source) const noexcept;

        [[nodiscard]] decltype(auto) calc(std::integral auto... dims) const;
        [[nodiscard]] decltype(auto) calc(IndexType index) const;
        [[nodiscard]] size_t toIndex1D(const IndexType& indices) const noexcept;
        [[nodiscard]] IndexType toIndexND(size_t index) const noexcept;
        void forND(std::invocable<T, IndexType> auto fn) const;

        void resize(const Tensor auto& x) { resize(x.getShape()); }
        decltype(auto) resize(std::integral auto... dims) { return Base::getDerived().resize(dims...); }
        decltype(auto) resize(IndexType shape) { return Base::getDerived().resize(shape); }

        [[nodiscard]] auto reals() const noexcept;
        [[nodiscard]] auto imags() const noexcept;
        [[nodiscard]] auto squaredNorms() const noexcept;
        [[nodiscard]] auto norms() const noexcept;
        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] constexpr int ndim() const noexcept { return NDim; }
        [[nodiscard]] size_t dim(int index) const noexcept;
        [[nodiscard]] IndexType getShape() const noexcept;
        [[nodiscard]] size_t getSize() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isForwardDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isReverseDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isLValueTensor() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isSparse() noexcept;
    protected:
        RValueTensor() = default;
        RValueTensor(const This&) = default;
        RValueTensor(This&&) noexcept = default;

        template<int GradOrder>
        auto grads_impl() const noexcept;
    };
}

namespace Physica {
    template<class T>
    class Traits<RValueTensor<T>> {
    public:
        using Derived = T;
    };
}

#include "RValueTensorImpl/RValueTensorImpl.h"
#include "RValueTensorImpl/TensorConvert.h"
#include "TensorExpr.h"
