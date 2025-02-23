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

namespace Physica {
    template<class Derived>
    template<Tensor T>
    Derived& LValueTensor<Derived>::operator=(const T& other) {
        if constexpr (std::is_same<Derived, T>::value)
            assert(this != &other && "[Error]: Self assign is likely a bug");
        resize(other.getDim());
        other.assign(Base::getDerived());
        return Base::getDerived();
    }

    template<class Derived>
    Derived& LValueTensor<Derived>::operator=(const ScalarType& s) {
        for (size_t i = 0; i < getDimX(); ++i)
            for (size_t j = 0; j < getDimY(); ++j)
                for (size_t k = 0; k < getDimZ(); ++k)
                    Base::getDerived()(i, j, k) = s;
        return Base::getDerived();
    }

    template<class Derived>
    inline LTensorBlock<Derived> LValueTensor<Derived>::block(Index3D from, Index3D count) {
        return {Base::getDerived(), from, count};
    }

    template<class Derived>
    inline const LTensorBlock<Derived> LValueTensor<Derived>::block(Index3D from, Index3D count) const {
        return {Base::getDerived(), from, count};
    }

    template<class Derived>
    template<RNG R>
    void LValueTensor<Derived>::random_uniform() {
        forIndexInTensor(getDim(), [this](Index3D index) {
            this->operator()(index) = ScalarType::template random_uniform<R>();
        });
    }

    template<class Derived>
    template<RNG R>
    void LValueTensor<Derived>::random_normal() {
        forIndexInTensor(getDim(), [this](Index3D index) {
            this->operator()(index) = ScalarType::template random_normal<R>();
        });
    }

    template<class Derived>
    void operator*=(LValueTensor<Derived>& grid, typename Derived::ScalarType factor) {
        for (size_t i = 0; i < grid.getDimX(); ++i)
            for (size_t j = 0; j < grid.getDimY(); ++j)
                for (size_t k = 0; k < grid.getDimZ(); ++k)
                    grid(i, j, k) *= factor;
    }
}
