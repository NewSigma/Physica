/*
 * Copyright 2023-2024 Weibo He.
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

#include <array>
#include <cassert>
#include "Physica/CRTPBase.h"
#include <H5Spublic.h>

namespace Physica {
    template<class Derived>
    class DataSpaceBase : public CRTPBase<DataSpaceBase<Derived>> {
        using This = DataSpaceBase<Derived>;
        using Base = CRTPBase<This>;
    public:
        constexpr static size_t Dim = Traits<Derived>::Dim;
        using SizeType = hsize_t;
        using SizeArray = std::array<hsize_t, Dim>;
    public:
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) { Base::getDerived().selectHyperslab(op, count, start); }
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return Base::getDerived().getDim(); }
        [[nodiscard]] size_t getSize(size_t dim) const noexcept { return Base::getDerived().getSize(dim); }
        [[nodiscard]] const SizeArray& getSelectedCount() const noexcept { return Base::getDerived().getSelectedCount(); }
        [[nodiscard]] const SizeArray& getSelectedStart() const noexcept { return Base::getDerived().getSelectedStart(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<DataSpaceBase<T>> {
    public:
        using Derived = T;
    };
}
