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

#include "../CompactTensor.h"

namespace Physica {
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
                stride *= Base::getDerived().dim(i);
            }
            return strides;
        }
    }
}
