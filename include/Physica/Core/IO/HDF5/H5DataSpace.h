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

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    class H5DataSpace : public H5::DataSpace {
        using Base = H5::DataSpace;
        using SizeArray = Utils::Array<hsize_t>;
    public:
        H5DataSpace(const H5::DataSpace& obj);
        H5DataSpace(const H5DataSpace&) = default;
        H5DataSpace(H5DataSpace&&) noexcept = delete;
        virtual ~H5DataSpace() = default;
        /* Operators */
        H5DataSpace& operator=(const H5DataSpace&) = default;
        H5DataSpace& operator=(H5DataSpace&&) noexcept = delete;
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start);
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return getSimpleExtentNdims(); }
        /* Static members */
        template<size_t Dim>
        [[nodiscard]] static inline H5DataSpace makeDataSpace(const Utils::Array<hsize_t, Dim>& dims);
    };

    template<size_t Dim>
    inline H5DataSpace H5DataSpace::makeDataSpace(const Utils::Array<hsize_t, Dim>& dims) {
        return H5DataSpace(H5::DataSpace(Dim == Utils::Dynamic ? dims.getLength() : Dim, dims.data()));
    }
}
