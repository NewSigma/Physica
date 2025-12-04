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

#include <cassert>
#include <sys/mman.h>
#include "Physica/Core/Exception/SystemException.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "HostAllocator.h"

namespace Physica {
    template<class T, size_t Align = Dynamic>
    class MMapAllocator : public HostAllocator<T, Align> {
        using This = MMapAllocator<T, Align>;
        using Base = HostAllocator<T, Align>;

        using Base::OverAlign;
        static_assert(!OverAlign, "[Error]: Not implemented");
    public:
        MMapAllocator() noexcept = default;
        MMapAllocator(const This&) noexcept = default;
        MMapAllocator(This&&) noexcept = delete;
        ~MMapAllocator() = default;
        /* Operators */
        This& operator=(const This&) noexcept = default;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard, gnu::returns_nonnull]] T* allocate(size_t n);
        void deallocate(T* p, size_t n);
        [[nodiscard]] T* reallocate(T* p, size_t new_size, size_t old_size);
        using Base::construct;
        using Base::destroy;
    };

    template<class T, size_t Align>
    T* MMapAllocator<T, Align>::allocate(size_t n) {
        assert(n > 0 && "[Error]: Allocate nothing");
        size_t size = Base::calcSize(n);
        constexpr int Prot = PROT_READ | PROT_WRITE;
        constexpr int Flag = MAP_PRIVATE;
        TempFile file(".MMAPXXXXXX");
        file.reserve(size);
        void* p = mmap(nullptr, size, Prot, Flag, file.getFd(), 0);
        if (p == MAP_FAILED)
            throw SystemException();
        return reinterpret_cast<T*>(p);
    }

    template<class T, size_t Align>
    void MMapAllocator<T, Align>::deallocate(T* p, size_t n) {
        if (munmap(p, Base::calcSize(n))) [[unlikely]]
            throw SystemException();
    }

    template<class T, size_t Align>
    T* MMapAllocator<T, Align>::reallocate(T* p, size_t new_size, size_t old_size) {
        assert(new_size > 0 && "[Error]: Allocate nothing");
        T* new_p = allocate(new_size);
        if (p == nullptr)
            return new_p;

        const size_t size = std::min(new_size, old_size);
        if constexpr (std::is_trivially_copyable<T>::value)
            memcpy(new_p, p, size * sizeof(T));
        else {
            for (size_t i = 0; i < size; ++i)
                construct(new_p + i, std::move(p[i]));
        }
        deallocate(p, old_size);
        return new_p;
    }
}

namespace std {
    template<class T, size_t Align>
    struct allocator_traits<Physica::MMapAllocator<T, Align>> : public allocator_traits<Physica::HostAllocator<T, Align>> {};
}
