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

#include <cstring>
#include "../CompactTensor.h"
#include "TensorBlock.h"
#include "TensorFiber.h"
#include "TensorSlice.h"

namespace Physica {
    template<class Derived>
    bool CompactTensor<Derived>::operator==(this const auto& self, const Tensor auto& other) noexcept {
        using U1 = std::remove_cvref_t<decltype(self)>;
        using U2 = std::remove_cvref_t<decltype(other)>;
        using T1 = U1::ScalarType;
        using T2 = U2::ScalarType;
        constexpr bool FastCompare = U1::isCompact() && U2::isCompact()
                                  && !U1::isDiffable() && !U2::isDiffable()
                                  && std::same_as<T1, T2>
                                  && std::is_trivially_copyable_v<T1>;
        if constexpr (FastCompare) {
            if (self.getShape() != other.getShape())
                return false;
            return std::memcmp(self.data(), other.data(), self.getSize() * sizeof(T)) == 0;
        }
        else
            return self.StridedTensor<Derived>::operator==(other);
    }

    template<class Derived>
    void CompactTensor<Derived>::assign(Tensor auto&& target) const noexcept {
        using Target = std::remove_cvref_t<decltype(target)>;
        constexpr bool FastCopy = Derived::isCompact() && Target::isCompact()
                               && !Derived::isDiffable() && !Target::isDiffable()
                               && std::same_as<typename Derived::ScalarType, typename Target::ScalarType>
                               && std::is_trivially_copyable_v<typename Derived::ScalarType>;
        const auto& self = Base::getDerived();
        target.assert_assign(self);
        if constexpr (FastCopy)
            std::memcpy(target.data(), self.data(), self.getSize() * sizeof(T));
        else
            Base::assign(target);
    }

    template<class Derived>
    auto CompactTensor<Derived>::data(this auto&& self) noexcept {
        return self.data_handle();
    }

    template<class Derived>
    auto CompactTensor<Derived>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto CompactTensor<Derived>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    constexpr auto CompactTensor<Derived>::getStrides() const noexcept {
        return getStridesImpl(Base::getDerived().getShape());
    }

    template<class Derived>
    __host__ __device__ constexpr auto CompactTensor<Derived>::getStridesImpl(const IndexType& shape) noexcept -> IndexType {
        constexpr bool KnownAtCompile = std::ranges::none_of(Derived::getStrideAtCompile(), [](size_t stride) consteval static {
            return stride == Dynamic;
        });

        if constexpr (KnownAtCompile)
            return Derived::getStrideAtCompile();
        else {
            IndexType strides = Derived::getStrideAtCompile();
            size_t stride = 1;
            for (int i = Base::NDim - 1; i >= 0; --i) {
                if (strides[i] == Dynamic)
                    strides[i] = stride;
                stride *= shape[i];
            }
            return strides;
        }
    }
}
