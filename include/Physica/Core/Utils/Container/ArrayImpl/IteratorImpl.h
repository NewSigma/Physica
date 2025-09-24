/*
 * Copyright 2025 Weibo He.
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

#include "Iterator.h"

namespace Physica {
    template<class Container>
    __host__ __device__ PtrIteratorF<Container>::PtrIteratorF(pointer p) noexcept : p(p) {}

    template<class Container>
    __host__ __device__ PtrIteratorF<Container>::PtrIteratorF(const This& ite) noexcept : p(ite.p) {}

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator+(difference_type n) const noexcept -> This {
        return This(p + n);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator+=(difference_type n) const noexcept -> This& {
        return *this = operator+(n);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator-(difference_type n) const noexcept -> This {
        return This(p - n);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator-=(difference_type n) const noexcept -> This& {
        return *this = operator-(n);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator-(const This& ite) const noexcept -> difference_type {
        return p - ite.p;
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorF<Container>::operator<(const This& ite) const noexcept {
        return p < ite.p;
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorF<Container>::operator>(const This& ite) const noexcept {
        return p > ite.p;
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorF<Container>::operator<=(const This& ite) const noexcept {
        return !(p > ite.p);
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorF<Container>::operator>=(const This& ite) const noexcept {
        return !(p < ite.p);
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorF<Container>::operator==(const This& ite) const noexcept {
        return p == ite.p;
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorF<Container>::operator!=(const This& ite) const noexcept {
        return p != ite.p;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator++() noexcept -> This& {
        ++p;
        return *this;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator++(int) noexcept -> This {
        return This(p++);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator--() noexcept -> This& {
        --p;
        return *this;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator--(int) noexcept -> This {
        return This(p--);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator*() const noexcept -> reference {
        return *p;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator->() const noexcept -> pointer {
        return p;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorF<Container>::operator[](difference_type n) const noexcept -> reference {
        return *operator+(n);
    }

    template<class Container>
    __host__ __device__ PtrIteratorR<Container>::PtrIteratorR(pointer p) noexcept : p(p) {}

    template<class Container>
    __host__ __device__ PtrIteratorR<Container>::PtrIteratorR(const This& ite) noexcept : p(ite.p) {}

    template<class Container>
    __host__ __device__ bool PtrIteratorR<Container>::operator==(const This& ite) const noexcept {
        return p == ite.p;
    }

    template<class Container>
    __host__ __device__ bool PtrIteratorR<Container>::operator!=(const This& ite) const noexcept {
        return p != ite.p;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorR<Container>::operator++() noexcept -> This&{
        --p; return *this;
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorR<Container>::operator++(int) noexcept -> This {
        return This(p--);
    }

    template<class Container>
    __host__ __device__ auto PtrIteratorR<Container>::operator*() const noexcept -> reference {
        return *p;
    }
}
