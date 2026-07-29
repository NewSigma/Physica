/*
 * Copyright 2023-2024 Weibo He.
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

#include <fstream>
#include <vector>
#include "H5Type.h"
#include "Mixin/Attributable.h"

namespace Physica {
    template<size_t Dim> class H5DataSpace;

    template<size_t Dim>
    class H5DataSet : public H5ID, public Attributable {
        using This = H5DataSet<Dim>;
    public:
        H5DataSet() = default;
        explicit H5DataSet(H5ID id_);
        H5DataSet(const This&) = default;
        H5DataSet(This&&) noexcept = default;
        ~H5DataSet() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<class MemSpace, class FileSpace>
        void read(void* buf, const H5Type& dtype, const MemSpace& mem_space, const FileSpace& file_space) const;
        void read(void* buf, const H5Type& dtype) const;

        template<class MemSpace, class FileSpace>
        void write(const void* buf, const H5Type& dtype, const MemSpace& mem_space, const FileSpace& file_space) const;
        void write(const void* buf, const H5Type& dtype) const;

        void readStr(char* buffer) const;
        void toFile(const char* path) const;
        /* Getters */
        [[nodiscard]] H5DataSpace<Dim> getDataSpace() const noexcept;
        [[nodiscard]] size_t getDim() const noexcept;
        [[nodiscard]] size_t getSize(size_t dim) const noexcept;
        [[nodiscard]] bool empty() const noexcept { return getSize(0) == 0; }
    };

    template<size_t Dim>
    H5DataSet<Dim>::H5DataSet(H5ID id_) : H5ID(std::move(id_)) {}

    template<size_t Dim>
    template<class MemSpace, class FileSpace>
    void H5DataSet<Dim>::read(void* buf, const H5Type& dtype, const MemSpace& mem_space, const FileSpace& file_space) const {
        H5Dread(getHID(), dtype.getHID(), mem_space.getHID(), file_space.getHID(), H5P_DEFAULT, buf);
    }

    template<size_t Dim>
    void H5DataSet<Dim>::read(void* buf, const H5Type& dtype) const {
        H5Dread(getHID(), dtype.getHID(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    }

    template<size_t Dim>
    template<class MemSpace, class FileSpace>
    void H5DataSet<Dim>::write(const void* buf, const H5Type& dtype, const MemSpace& mem_space, const FileSpace& file_space) const {
        H5Dwrite(getHID(), dtype.getHID(), mem_space.getHID(), file_space.getHID(), H5P_DEFAULT, buf);
    }

    template<size_t Dim>
    void H5DataSet<Dim>::write(const void* buf, const H5Type& dtype) const {
        H5Dwrite(getHID(), dtype.getHID(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    }

    template<size_t Dim>
    void H5DataSet<Dim>::readStr(char* buffer) const {
        H5Dread(getHID(), H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        buffer[getSize(0)] = '\0';
    }

    template<size_t Dim>
    void H5DataSet<Dim>::toFile(const char* path) const {
        std::ofstream fout(path);
        const auto size = getDataSpace().getSize(0);
        auto buffer = std::vector<char>(size);
        fout.write(buffer.data(), size);
    }

    template<size_t Dim>
    H5DataSpace<Dim> H5DataSet<Dim>::getDataSpace() const noexcept {
        return H5DataSpace<Dim>(H5ID(H5Dget_space(getHID())));
    }

    template<size_t Dim>
    size_t H5DataSet<Dim>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return getDataSpace().getDim();
        else
            return Dim;
    }

    template<size_t Dim>
    size_t H5DataSet<Dim>::getSize(size_t dim) const noexcept {
        return getDataSpace().getSize(dim);
    }
}
