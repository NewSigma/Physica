/*
 * Copyright 2025-2026 Weibo He.
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

#include "../RValueTensor.h"
#include "Physica/Core/Utils/Container/ArrayND.h"

namespace Physica {
    template<class Derived, Scalar ScalarT>
    void RValueTensor<Derived, ScalarT>::assign(Tensor auto& x) const {
        x.assert_assign(Base::getDerived());
        if constexpr (!isDiffable() && x.isDiffable()) {
            Base::getDerived().assign(x.values());
            x.zero_grad();
        }
        else {
            size_t size = Base::getDerived().getSize();
            for (size_t i = 0; i < size; ++i) {
                const auto indices = toIndexND(i);
                x[indices] = calc(indices);
            }
        }
    }

    template<class Derived, Scalar ScalarT>
    void RValueTensor<Derived, ScalarT>::assert_assign(const Tensor auto& source) const noexcept {
        if constexpr (std::is_same<Derived, std::remove_cvref_t<decltype(source)>>::value)
            assert(this != &source && "[Error]: Self assign is likely a bug");
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueTensor<Derived, ScalarT>::calc(std::integral auto... dims) const {
        return Base::getDerived().calc(dims...);
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueTensor<Derived, ScalarT>::calc(IndexType index) const {
        return Base::getDerived().calc(index);
    }

    template<class Derived, Scalar ScalarT>
    size_t RValueTensor<Derived, ScalarT>::toIndex1D(const IndexType& indices) const noexcept {
        return IndexType::toIndex1D(getShape(), indices);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueTensor<Derived, ScalarT>::toIndexND(size_t index) const noexcept -> IndexType {
        return IndexType::toIndexND(getShape(), index);
    }

    template<class Derived, Scalar ScalarT>
    void RValueTensor<Derived, ScalarT>::forND(std::invocable<T, IndexType> auto fn) const {
        Physica::forND(getShape(), [this, fn](const IndexType& index) {
            fn(calc(index), index);
        });
    }

    template<class Derived, Scalar ScalarT>
    auto RValueTensor<Derived, ScalarT>::reals() const noexcept {
        return RealTensor<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueTensor<Derived, ScalarT>::imags() const noexcept {
        return ImagTensor<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueTensor<Derived, ScalarT>::squaredNorms() const noexcept {
        return SquaredNormTensor<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueTensor<Derived, ScalarT>::norms() const noexcept {
        return NormTensor<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueTensor<Derived, ScalarT>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return ValueVector<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    template<int GradOrder>
    auto RValueTensor<Derived, ScalarT>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        return GradTensor<Derived, GradOrder>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    size_t RValueTensor<Derived, ScalarT>::dim(int index) const noexcept {
        return Base::getDerived().dim(index);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueTensor<Derived, ScalarT>::getShape() const noexcept -> IndexType {
        return Base::getDerived().getShape();
    }

    template<class Derived, Scalar ScalarT>
    size_t RValueTensor<Derived, ScalarT>::getSize() const noexcept {
        size_t size = dim(0);
        for (int i = 1; i < NDim; ++i)
            size *= dim(i);
        return size;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isDiffable() noexcept {
        return Diffable<T>;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isLValueTensor() noexcept {
        return false;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isCompact() noexcept {
        return requires{ std::declval<Derived>().data(); };
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueTensor<Derived, ScalarT>::isSparse() noexcept {
        return requires{ std::declval<Derived>().getNumNonzero(); };
    }

    template<class Derived, Scalar ScalarT>
    template<IndexVar... Ts>
    __host__ __device__ consteval int RValueTensor<Derived, ScalarT>::calcFiberDim() noexcept {
        constexpr IndexVarInfo<Ts...> info{};
        static_assert(info.getNumAnonymous() == 1, "[Error]: Fiber requires 1 anonymous var");
        return std::ranges::distance(info.isAnonymous.begin(), std::ranges::find(info.isAnonymous, true));
    }

    template<class Derived, Scalar ScalarT>
    template<IndexVar... Ts>
    __host__ __device__ consteval Array<int, 2> RValueTensor<Derived, ScalarT>::calcSliceDim() noexcept {
        constexpr IndexVarInfo<Ts...> info{};
        static_assert(info.getNumAnonymous() == 2, "[Error]: Slice requires 2 anonymous var");
        Array<int, 2> result{};
        int count = 0;
        for (int i = 0; i < info.size() && count < 2; ++i) {
            if (info.isAnonymous[i])
                result[count++] = i;
        }
        return result;
    }
}
