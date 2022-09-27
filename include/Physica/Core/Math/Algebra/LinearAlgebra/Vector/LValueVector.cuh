/*
 * Copyright 2022 WeiBo He.
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

#include "LValueVector.h"

namespace Physica::Core {
    template<class Derived>
    class device_obj<LValueVector<Derived>> : public device_obj<RValueVector<Derived>> {
    public:
        using Base = device_obj<RValueVector<Derived>>;
        using typename Base::ScalarType;
    public:
        /* Operators */
        device_obj<LValueVector<Derived>>& operator=(const device_obj<LValueVector<Derived>>& obj);
        device_obj<LValueVector<Derived>>& operator=(device_obj<LValueVector<Derived>>&& obj) noexcept;
        template<class T>
        device_obj<Derived>& operator=(const device_obj<RValueVector<T>>& v);
        [[nodiscard]] __device__ ScalarType& operator[](size_t index) { return Base::getDerived()[index]; }
        [[nodiscard]] __device__ const ScalarType& operator[](size_t index) const { return Base::getDerived()[index]; }
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return (*this)[index]; }
    };

    template<class Derived>
    device_obj<LValueVector<Derived>>& device_obj<LValueVector<Derived>>::operator=(const device_obj<LValueVector<Derived>>& obj) {
        return Base::getDerived() = obj.getDerived();
    }
    template<class Derived>
    device_obj<LValueVector<Derived>>& device_obj<LValueVector<Derived>>::operator=(device_obj<LValueVector<Derived>>&& obj) noexcept {
        return Base::getDerived() = std::move(obj.getDerived());
    }

    template<class Derived>
    template<class T>
    device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const device_obj<RValueVector<T>>& v) {
        Base::getDerived().resize(v.getLength());
        v.assignTo(*this);
        return *this;
    }
}
