/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <array>
#include <cstdlib>
#include <unistd.h>
#include "Physica/Core/Exception/Syscall.h"

namespace Physica {
    template<size_t N>
    class TempFile {
        using This = TempFile<N>;

        std::array<char, N> name;
        int fd;
    public:
        TempFile(const char (&nameFmt)[N]);
        TempFile(const This&) = delete;
        TempFile(This&& obj) noexcept;
        ~TempFile();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void reserve(int64_t size);
        void release() noexcept { fd = -1; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const char* getName() const noexcept { return name.data(); }
        [[nodiscard]] int getFd() const noexcept { return fd; }
    };

    template<size_t N>
    TempFile<N>::TempFile(const char (&nameFmt)[N])
            : name(std::to_array(nameFmt))
            , fd(mkstemp(name.data())) {
        check_syscall(fd == -1);
    }

    template<size_t N>
    TempFile<N>::TempFile(This&& obj) noexcept : fd(obj.fd) {
        for (size_t i = 0; i < N; ++i)
            name[i] = obj.name[i];
        obj.fd = -1;
    }

    template<size_t N>
    TempFile<N>::~TempFile() {
        if (fd != -1) {
            close(fd);
            unlink(name.data());
            fd = -1;
        }
    }

    template<size_t N>
    void TempFile<N>::reserve(int64_t size) {
        check_syscall(ftruncate64(fd, size) == -1);
    }

    template<size_t N>
    void TempFile<N>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(name, obj.name);
        std::swap(fd, obj.fd);
    }
}
