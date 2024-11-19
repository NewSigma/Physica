/*
 * Copyright 2023 Weibo He.
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

namespace Physica::Core {
    template<class Derived>
    template<Grid T>
    Derived& LValueGrid<Derived>::operator=(const T& other) {
        if constexpr (std::is_same<Derived, T>::value)
            assert(this != &other && "[Error]: Self assign is likely a bug");
        resize(other.getDim());
        other.assignTo(*this);
        return Base::getDerived();
    }

    template<class Derived>
    Derived& LValueGrid<Derived>::operator=(const ScalarType& s) {
        for (size_t i = 0; i < getDimX(); ++i)
            for (size_t j = 0; j < getDimY(); ++j)
                for (size_t k = 0; k < getDimZ(); ++k)
                    Base::getDerived()(i, j, k) = s;
        return Base::getDerived();
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::leftFrontBottomCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {0, 0, 0}, cornerIndex};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::leftFrontBottomCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {0, 0, 0}, cornerIndex};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::rightFrontBottomCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {cornerIndex[0], 0, 0}, {getDimX() - cornerIndex[0], cornerIndex[1], cornerIndex[2]}};
    }
    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::rightFrontBottomCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {cornerIndex[0], 0, 0}, {getDimX() - cornerIndex[0], cornerIndex[1], cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::leftBackBottomCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {0, cornerIndex[1], 0}, {cornerIndex[0], getDimY() - cornerIndex[1], cornerIndex[2]}};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::leftBackBottomCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {0, cornerIndex[1], 0}, {cornerIndex[0], getDimY() - cornerIndex[1], cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::rightBackBottomCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {cornerIndex[0], cornerIndex[1], 0}, {getDimX() - cornerIndex[0], getDimY() - cornerIndex[1], cornerIndex[2]}};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::rightBackBottomCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {cornerIndex[0], cornerIndex[1], 0}, {getDimX() - cornerIndex[0], getDimY() - cornerIndex[1], cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::leftFrontTopCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {0, 0, cornerIndex[2]}, {cornerIndex[0], cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::leftFrontTopCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {0, 0, cornerIndex[2]}, {cornerIndex[0], cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::rightFrontTopCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {cornerIndex[0], 0, cornerIndex[2]}, {getDimX() - cornerIndex[0], cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }
    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::rightFrontTopCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {cornerIndex[0], 0, cornerIndex[2]}, {getDimX() - cornerIndex[0], cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::leftBackTopCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {0, cornerIndex[1], cornerIndex[2]}, {cornerIndex[0], getDimY() - cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::leftBackTopCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {0, cornerIndex[1], cornerIndex[2]}, {cornerIndex[0], getDimY() - cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::rightBackTopCorner(Index3D cornerIndex) {
        return {Base::getDerived(), {cornerIndex[0], cornerIndex[1], cornerIndex[2]}, {getDimX() - cornerIndex[0], getDimY() - cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::rightBackTopCorner(Index3D cornerIndex) const {
        return {Base::getConstCastDerived(), {cornerIndex[0], cornerIndex[1], cornerIndex[2]}, {getDimX() - cornerIndex[0], getDimY() - cornerIndex[1], getDimZ() - cornerIndex[2]}};
    }

    template<class Derived>
    inline LGridBlock<Derived> LValueGrid<Derived>::block(Index3D from, Index3D count) {
        return {Base::getDerived(), from, count};
    }

    template<class Derived>
    inline const LGridBlock<Derived> LValueGrid<Derived>::block(Index3D from, Index3D count) const {
        return {Base::getDerived(), from, count};
    }

    template<class Derived>
    template<class RandomType>
    void LValueGrid<Derived>::random_uniform(RandomType& gen) {
        forIndexInGrid(getDim(), [this, &gen](Index3D index) {
            this->operator()(index) = ScalarType::random_uniform(gen);
        });
    }

    template<class Derived>
    template<class RandomType>
    void LValueGrid<Derived>::random_normal(RandomType& gen) {
        forIndexInGrid(getDim(), [this, &gen](Index3D index) {
            this->operator()(index) = ScalarType::random_normal(gen);
        });
    }

    template<class Derived>
    void operator*=(LValueGrid<Derived>& grid, typename Derived::ScalarType factor) {
        for (size_t i = 0; i < grid.getDimX(); ++i)
            for (size_t j = 0; j < grid.getDimY(); ++j)
                for (size_t k = 0; k < grid.getDimZ(); ++k)
                    grid(i, j, k) *= factor;
    }
}
