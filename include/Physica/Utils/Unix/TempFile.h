/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Utils {
    template<size_t N>
    class TempFile {
        char name[N];
        int fd;
    public:
        TempFile(const char (&name_template)[N]);
        TempFile(const TempFile&) = delete;
        TempFile(TempFile&& obj) noexcept;
        ~TempFile();
        /* Operators */
        TempFile& operator=(TempFile obj) noexcept;
        /* Operations */
        void swap(TempFile& obj) noexcept;
        /* Getters */
        [[nodiscard]] const char* getName() const noexcept { return name; }
        [[nodiscard]] int getFd() const noexcept { return fd; }
    };

    template<size_t N>
    TempFile<N>::TempFile(const char (&name_template)[N]) {
        for (size_t i = 0; i < N; ++i)
            name[i] = name_template[i];
        fd = mkstemp(name);
    }

    template<size_t N>
    TempFile<N>::TempFile(TempFile&& obj) noexcept : fd(obj.fd) {
        for (size_t i = 0; i < N; ++i)
            name[i] = obj.name[i];
        obj.fd = -1;
    }

    template<size_t N>
    TempFile<N>::~TempFile() {
        if (fd != -1) {
            close(fd);
            unlink(name);
            fd = -1;
        }
    }

    template<size_t N>
    TempFile<N>& TempFile<N>::operator=(TempFile<N> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<size_t N>
    void TempFile<N>::swap(TempFile& obj) noexcept {
        std::swap(name, obj.name);
        std::swap(fd, obj.fd);
    }
}
