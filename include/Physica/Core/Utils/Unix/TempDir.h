/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstdlib>
#include <filesystem>
#include <unistd.h>
#include "Physica/Core/Exception/Syscall.h"

namespace Physica {
    template<size_t N>
    class TempDir {
        using This = TempDir<N>;

        char name[N];
        const char* pName;
    public:
        TempDir(const char (&name_template)[N]);
        TempDir(const TempDir&) = delete;
        TempDir(TempDir&& obj) noexcept;
        ~TempDir();
        /* Operators */
        TempDir& operator=(TempDir obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const char* getName() const noexcept { return name; }
    };

    template<size_t N>
    TempDir<N>::TempDir(const char (&name_template)[N]) {
        for (size_t i = 0; i < N; ++i)
            name[i] = name_template[i];
        pName = mkdtemp(name);
        check_syscall(pName == nullptr);
    }

    template<size_t N>
    TempDir<N>::TempDir(TempDir&& obj) noexcept : pName(obj.pName) {
        for (size_t i = 0; i < N; ++i)
            name[i] = obj.name[i];
        obj.pName = nullptr;
    }

    template<size_t N>
    TempDir<N>::~TempDir() {
        if (pName != nullptr) {
            std::error_code error;
            std::filesystem::remove_all(pName, error);
            pName = nullptr;
        }
    }

    template<size_t N>
    void TempDir<N>::swap(TempDir& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(name, obj.name);
        std::swap(pName, obj.pName);
    }
}
