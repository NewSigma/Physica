/*
 * Copyright 2024-2025 Weibo He.
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

#include "../ContinuousVector.h"

namespace Physica {
    template<class T>
    class ValueVector<ContinuousVector<T>> : public ContinuousVector<ValueVector<ContinuousVector<T>>> {
        using This = ValueVector<ContinuousVector<T>>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        const T& v;
    public:
        ValueVector(const T& v_) : v(v_) {}
        ValueVector(const This&) = default;
        ValueVector(This&&) noexcept = default;
        ~ValueVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;

        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == v.getLength()); }
        /* Getters */
        [[nodiscard]] PtrTy data_ptr(size_t i) noexcept { return (*v.data_ptr(i)).value_ptr(); }
        [[nodiscard]] ConstPtrTy data_ptr(size_t i) const noexcept { return const_cast<This&>(*this).data_ptr(i); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class T, int GradOrder>
    class GradVector<ContinuousVector<T>, GradOrder> : public ContinuousVector<GradVector<ContinuousVector<T>, GradOrder>> {
        using This = GradVector<ContinuousVector<T>, GradOrder>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        const T& v;
    public:
        GradVector(const T& v_) : v(v_) {}
        GradVector(const This&) = default;
        GradVector(This&&) noexcept = default;
        ~GradVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;

        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == v.getLength()); }
        /* Getters */
        [[nodiscard]] PtrTy data_ptr(size_t i) noexcept { return v.data_ptr(i).template grad_ptr<GradOrder>(); }
        [[nodiscard]] ConstPtrTy data_ptr(size_t i) const noexcept { return const_cast<This&>(*this).data_ptr(i); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };
}
