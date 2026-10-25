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

#include <algorithm>
#include <utility>
#include "../StridedTensor.h"
#include "TensorFiber.h"

namespace Physica {
    template<class Derived>
    bool StridedTensor<Derived>::operator==(this const auto& self, const Tensor auto& other) noexcept {
        if constexpr (isMinorCompact<decltype(self), decltype(other)>()) {
            if (self.getShape() != other.getShape())
                return false;
            const size_t size = self.getSize();
            const size_t stride = std::min(getCompactMinorSize(self), getCompactMinorSize(other));
            for (size_t i = 0; i < size; i += stride) {
                const auto indices = self.toIndexND(i);
                if (std::memcmp(self.data_ptr(indices), other.data_ptr(indices), stride * sizeof(T)) != 0)
                    return false;
            }
            return true;
        }
        return self.LValueTensor<Derived>::operator==(other);
    }

    template<class Derived>
    void StridedTensor<Derived>::assign(Tensor auto&& target) const noexcept {
        using Target = std::remove_cvref_t<decltype(target)>;
        if constexpr (isMinorCompact<Derived, Target>()) {
            const auto& source = Base::getDerived();
            const int minorDim = std::max(getCompactMinorDim(source), getCompactMinorDim(target));
            const size_t size = source.getSize();
            const size_t stride = std::min(getCompactMinorSize(source), getCompactMinorSize(target));
            for (size_t i = 0; i < size; i += stride) {
                const auto indices = source.toIndexND(i);
                auto src = makeCompactBlock(source, indices, minorDim);
                auto tgt = makeCompactBlock(target, indices, minorDim);
                src.flatten().assign(tgt.flatten());
            }
        }
        else
            Base::assign(target);
    }

    template<class Derived>
    void StridedTensor<Derived>::zeros() noexcept {
        if constexpr (isMinorCompact<Derived, Derived>()) {
            auto& self = Base::getDerived();
            const int minorDim = getCompactMinorDim(self);
            const size_t size = self.getSize();
            const size_t stride = getCompactMinorSize(self);
            for (size_t i = 0; i < size; i += stride) {
                const auto indices = self.toIndexND(i);
                makeCompactBlock(self, indices, minorDim).flatten().zeros();
            }
        }
        else
            Base::zeros();
    }

    template<class Derived>
    auto StridedTensor<Derived>::getStrides() const noexcept {
        return Base::getDerived().getStrides();
    }

    template<class Derived>
    size_t StridedTensor<Derived>::getStride(size_t dim) const noexcept {
        return getStrides()[dim];
    }

    template<class Derived>
    auto StridedTensor<Derived>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedTensor<Derived>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedTensor<Derived>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        const IndexType strides = self.getStrides();
        size_t offset = 0;
        for (int i = 0; i < NDim; ++i) {
            assert(index[i] < self.dim(i));
            offset += index[i] * strides[i];
        }
        return self.data_handle() + offset;
    }

    template<class Derived>
    auto StridedTensor<Derived>::data_ptr(this auto&& self, std::same_as<size_t> auto... dims) noexcept {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        return self.data_ptr(IndexType({dims...}));
    }

    template<class Derived>
    __host__ __device__ consteval auto StridedTensor<Derived>::getStrideAtCompile() noexcept -> IndexType {
        return IndexType::zeros(NDim);
    }

    template<class Derived>
    __host__ __device__ consteval auto StridedTensor<Derived>::getStrideAtCompile(size_t dim) noexcept -> size_t {
        return Derived::getStrideAtCompile()[dim];
    }

    template<class Derived>
    template<class Self, class Other>
    __host__ __device__ consteval bool StridedTensor<Derived>::isMinorCompact() noexcept {
        using U1 = std::remove_cvref_t<Self>;
        using U2 = std::remove_cvref_t<Other>;
        using T1 = U1::ScalarType;
        using T2 = U2::ScalarType;
        if constexpr (U1::isStrided() && U2::isStrided())
            return U1::getStrideAtCompile().back() == 1
                && U2::getStrideAtCompile().back() == 1
                && !U1::isDiffable() && !U2::isDiffable()
                && std::same_as<T1, T2>
                && std::is_trivially_copyable_v<T1>;
        else
            return false;
    }

    template<class Derived>
    __host__ __device__ constexpr int StridedTensor<Derived>::getCompactMinorDim(const auto& tensor) noexcept {
        const auto shape = tensor.getShape();
        const auto strides = tensor.getStrides();
        int dim = tensor.ndim();
        size_t result = 1;
        for (int i = tensor.ndim() - 1; i >= 0; --i) {
            if (strides[i] != result)
                break;
            result *= shape[i];
            dim = i;
        }
        return dim;
    }

    template<class Derived>
    __host__ __device__ constexpr size_t StridedTensor<Derived>::getCompactMinorSize(const auto& tensor) noexcept {
        const auto shape = tensor.getShape();
        size_t result = 1;
        for (int i = getCompactMinorDim(tensor); i < tensor.ndim(); ++i)
            result *= shape[i];
        return result;
    }

    template<class Derived>
    __host__ __device__ constexpr auto StridedTensor<Derived>::makeCompactBlock(auto&& tensor, const IndexType& index, int minorDim) noexcept {
        IndexType count = tensor.getShape();
        for (int i = 0; i < minorDim; ++i)
            count[i] = 1;
        return tensor.block(index, std::move(count));
    }
}
