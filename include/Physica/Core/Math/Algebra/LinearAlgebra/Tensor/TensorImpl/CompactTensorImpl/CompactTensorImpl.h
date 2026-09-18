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
    auto CompactTensor<Derived>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        return self.data_handle() + self.toIndex1D(index);
    }

    template<class Derived>
    auto CompactTensor<Derived>::data_ptr(this auto&& self, std::integral auto... dims) noexcept {
        static_assert(sizeof...(dims) == Base::NDim, "[Error]: NDim is not consistent");
        return self.data_ptr(IndexType({static_cast<size_t>(dims)...}));
    }
}
