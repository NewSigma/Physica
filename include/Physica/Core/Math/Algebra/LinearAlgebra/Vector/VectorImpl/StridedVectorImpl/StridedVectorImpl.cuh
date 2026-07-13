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

#include "../StridedVector.cuh"

namespace Physica {
    template<class Derived>
    auto device_obj<StridedVector<Derived>>::operator=(Scalar auto x) -> device_obj<Derived>& {
        if constexpr (Base::isCompact() && Base::getSizeAtCompile() == Dynamic) {
            if (x.isZero())
                zeros();
        }
        return Base::operator=(x);
    }

    template<class Derived>
    void device_obj<StridedVector<Derived>>::zeros() {
        if constexpr (Base::isCompact()) {
            if constexpr (Diffable<T>) {
                Base::getDerived().values().zeros();
                Base::getDerived().zero_grad();
            }
            else
                check(cudaMemsetAsync(data_handle(), 0, Base::getLength() * sizeof(T)));
        }
        else
            Base::zeros();
    }

    template<class Derived>
    template<RNG R>
    void device_obj<StridedVector<Derived>>::random_uniform() {
        if constexpr (Base::isCompact() && R::cuRAND_Ready) {
            auto& rng = R::getInstance();
            check(curandSetStream(rng, CUDAContext::getInstance()));
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1);
            if constexpr (T::Prec == Float32)
                check(curandGenerateUniform(rng, (float*)data_handle(), length));
            else if constexpr (T::Prec == Float64)
                check(curandGenerateUniformDouble(rng, (double*)data_handle(), length));
            else
                Base::template random_uniform<R>();
        }
        else
            Base::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void device_obj<StridedVector<Derived>>::random_normal() {
        if constexpr (Base::isCompact() && R::cuRAND_Ready) {
            auto& rng = R::getInstance();
            check(curandSetStream(rng, CUDAContext::getInstance()));
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1);
            if constexpr (T::Prec == Float32)
                check(curandGenerateNormal(rng, (float*)data_handle(), length, 0, 1));
            else if constexpr (T::Prec == Float64)
                check(curandGenerateNormalDouble(rng, (double*)data_handle(), length, 0, 1));
            else
                Base::template random_normal<R>();
        }
        else
            Base::template random_normal<R>();
    }
}
