/*
 * Copyright 2024 Weibo He.
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
    template<Matrix M>
    __host__ __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const device_obj<M>& m) {
        static_assert(RowAtCompile == Dynamic || M::RowAtCompile == Dynamic || RowAtCompile == M::RowAtCompile, "[Error]: Incompatible row number");
        static_assert(ColAtCompile == Dynamic || M::ColAtCompile == Dynamic || ColAtCompile == M::ColAtCompile, "[Error]: Incompatible col number");
        auto& target = Base::getDerived();
        target.resize(m.getRow(), m.getCol());
        m.assign(target);
        return target;
    }

    template<class Derived>
    __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const ScalarType& x) {
        for (size_t i = 0; i < Base::getMaxMajor(); ++i)
            for (size_t j = 0; j < Base::getMaxMinor(); ++j)
                refFromMajorMinor(i, j) = ScalarType(x);
        return Base::getDerived();
    }

    template<class Derived>
    __device__ inline device_obj<LValueMatrix<Derived>>::ScalarType&
    device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::colFromMajorMinor<Derived>(major, minor);
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getCol());
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    __device__ inline const device_obj<LValueMatrix<Derived>>::ScalarType&
    device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) const {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::colFromMajorMinor<Derived>(major, minor);
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getCol());
        return Base::getDerived()(r, c);
    }
}
