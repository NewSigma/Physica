/*
 * Copyright 2024 WeiBo He.
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

#include <cstdlib>
#include <unistd.h>
#include "UnixHelper.h"
#include "Physica/Core/Exception/SyscallException.h"

namespace Physica::Utils {
    template<size_t N>
    class TempDir {
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
        void swap(TempDir& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const char* getName() const noexcept { return name; }
    };

    template<size_t N>
    TempDir<N>::TempDir(const char (&name_template)[N]) {
        for (size_t i = 0; i < N; ++i)
            name[i] = name_template[i];
        pName = mkdtemp(name);
        if (pName == nullptr) [[unlikely]]
            throw Core::SyscallException();
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
            forceRemoveDir(pName);
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
