/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cassert>
#include <sys/mman.h>
#include "Physica/Core/Exception/Syscall.h"
#include "Physica/Core/Utils/Builtin.h"
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
        MMapAllocator(This&&) noexcept = default;
        ~MMapAllocator() = default;
        /* Operators */
        This& operator=(const This&) noexcept = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::returns_nonnull]] T* allocate(size_t n);
        void deallocate(T* p, size_t n);
        [[nodiscard, gnu::returns_nonnull]] T* reallocate(T* p, size_t new_size, size_t old_size);
    protected:
        [[nodiscard]] static size_t calcSize(size_t n) noexcept;
    };

    template<class T, size_t Align>
    T* MMapAllocator<T, Align>::allocate(size_t n) {
        size_t size = calcSize(n);
        constexpr int Prot = PROT_READ | PROT_WRITE;
        constexpr int Flag = MAP_PRIVATE;
        TempFile file(".MMAPXXXXXX");
        file.reserve(size);
        void* p = mmap(nullptr, size, Prot, Flag, file.getFd(), 0);
        check_syscall(p == MAP_FAILED);
        return reinterpret_cast<T*>(p);
    }

    template<class T, size_t Align>
    void MMapAllocator<T, Align>::deallocate(T* p, size_t n) {
        check_syscall(munmap(p, calcSize(n)) == -1);
    }

    template<class T, size_t Align>
    T* MMapAllocator<T, Align>::reallocate(T* p, size_t new_size, size_t old_size) {
        assert(new_size > 0 && "[Error]: Reject bad pattern");
        assert(p != nullptr || old_size == 0);
        T* new_p = allocate(new_size);
        memcpy((void*)new_p, p, std::min(new_size, old_size) * sizeof(T));
        deallocate(p, old_size);
        return new_p;
    }

    template<class T, size_t Align>
    size_t MMapAllocator<T, Align>::calcSize(size_t n) noexcept {
        size_t size = ckd_mul(n, sizeof(T));
        if constexpr (OverAlign)
            size = (size + Align - 1) & ~(Align - 1);
        return size;
    }
}

namespace std {
    template<class T, size_t Align>
    struct allocator_traits<Physica::MMapAllocator<T, Align>> : public allocator_traits<Physica::HostAllocator<T, Align>> {};
}
