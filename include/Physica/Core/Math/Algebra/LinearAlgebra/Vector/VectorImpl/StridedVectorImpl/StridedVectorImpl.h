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

#include "../StridedVector.h"
#include "Physica/Core/Exception/MKL/VSL.h"

namespace Physica {
    template<class Derived>
    auto StridedVector<Derived>::operator=(Scalar auto x) noexcept -> Derived& {
        if constexpr (isCompact() && Base::getSizeAtCompile() == Dynamic) {
            if (x.isZero())
                zeros();
        }
        return Base::operator=(x);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void StridedVector<Derived>::assign(Vector auto&& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        constexpr bool SameScalar = std::same_as<T, typename V::ScalarType>;
        constexpr bool Copyable = std::is_trivially_copyable<T>::value;
        if constexpr (isCompact() && v.isCompact() && SameScalar && Copyable) {
            if constexpr (Base::isDiffable()) {
                auto& x = Base::getDerived();
                x.values().template assign<P>(v.values());
                x.grads().template assign<P>(v.grads());
                return;
            }
            else {
                constexpr size_t Size = std::max(Base::getSizeAtCompile(), v.getSizeAtCompile());
                if constexpr (Size != Dynamic)
                    memcpy(v.data(), Base::getDerived().data(), Size * sizeof(T)); // Static memcpy
                else {
                    constexpr size_t Critical = 1024 * sizeof(float64); // Based on benchmark
                    constexpr bool MaybeBenefit = Size == Dynamic || (Size * sizeof(T) > Critical);
                    size_t size = Base::getLength() * sizeof(T);
                    if (MaybeBenefit && (size > Critical))
                        memcpy(v.data(), Base::getDerived().data(), size); // Reuse size, because v.read(data()) shows IR regression
                    else
                        Base::template assign_base<P>(v);
                }
            }
        }
        else
            Base::template assign_base<P>(v);
    }

    template<class Derived>
    auto StridedVector<Derived>::norm1() const noexcept -> CoDiff<Tr> {
        constexpr size_t Size = Base::getSizeAtCompile();
        constexpr bool SmallVector = 0 < Size && Size <= 128;
        if constexpr (HasMKL() && Internal::EnableLAPACK<Derived>::value && !SmallVector) {
            bool isSmallVector = Base::getLength() <= 128;
            return isSmallVector ? norm1_base() : norm1_mkl();
        }
        else
            return norm1_base();
    }

    template<class Derived>
    auto StridedVector<Derived>::norm1_base() const noexcept -> CoDiff<Tr> {
        return Base::norm1();
    }

    template<class Derived>
    auto StridedVector<Derived>::norm2() const noexcept -> CoDiff<Tr> {
        return norm2_base();
    }

    template<class Derived>
    auto StridedVector<Derived>::norm2_base() const noexcept -> CoDiff<Tr> {
        return Base::norm2();
    }
    /**
     * Prefer zeros() over simply assigning zero for better performance.
     */
    template<class Derived>
    void StridedVector<Derived>::zeros() noexcept {
        if constexpr (isCompact()) {
            if constexpr (Diffable<T>) {
                Base::getDerived().values().zeros();
                Base::getDerived().zero_grad();
            }
            else if constexpr (T::Prec != FloatMP)
                std::memset(data_handle(), 0, Base::getLength() * sizeof(T));
            else
                Base::zeros();
        }
        else
            Base::zeros();
    }

    template<class Derived>
    template<RNG R>
    void StridedVector<Derived>::random_uniform() {
        if constexpr (isCompact() && R::MKL_Ready) {
            [[maybe_unused]] auto& gen = R::getInstance();
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1) * (Base::isForwardDiff() ? 2 : 1);
            if constexpr (T::Prec == Float32)
                check_vsl(vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (float*)data_handle(), 0, 1));
            else if constexpr (T::Prec == Float64)
                check_vsl(vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (double*)data_handle(), 0, 1));
            else
                Base::template random_uniform<R>();
        }
        else
            Base::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void StridedVector<Derived>::random_normal() {
        if constexpr (isCompact() && R::MKL_Ready && !Base::isForwardDiff()) {
            [[maybe_unused]] auto& gen = R::getInstance();
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1);
            if constexpr (T::Prec == Float32)
                check_vsl(vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (float*)data_handle(), 0, 1));
            else if constexpr (T::Prec == Float64)
                check_vsl(vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (double*)data_handle(), 0, 1));
            else
                Base::template random_normal<R>();
        }
        else
            Base::template random_normal<R>();
    }
}
