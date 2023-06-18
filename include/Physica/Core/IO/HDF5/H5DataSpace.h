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
    template<size_t Dim>
    class H5DataSpace : public H5::DataSpace {
        using Base = H5::DataSpace;
        using SizeArray = Utils::Array<hsize_t, Dim>;
        constexpr static auto Dynamic = Physica::Utils::Dynamic;
    public:
        H5DataSpace(const H5::DataSpace& obj);
        H5DataSpace(const H5DataSpace&) = default;
        H5DataSpace(H5DataSpace&&) noexcept = delete;
        virtual ~H5DataSpace() = default;
        /* Operators */
        inline H5DataSpace& operator=(const H5DataSpace& obj);
        H5DataSpace& operator=(H5DataSpace&&) noexcept = delete;
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start);
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept;
        [[nodiscard]] size_t getSize(size_t dim) const noexcept;
        /* Static members */
        [[nodiscard]] static inline H5DataSpace makeDataSpace(const SizeArray& dims);
        [[nodiscard]] static inline H5DataSpace makeDataSpace(const SizeArray& dims, const SizeArray& maxdims);
    };

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(const H5::DataSpace& obj) : Base(obj) {}

    template<size_t Dim>
    void H5DataSpace<Dim>::selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) {
        assert(count.getLength() == getDim());
        assert(start.getLength() == getDim());
        Base::selectHyperslab(op, count.data(), start.data());
    }

    template<size_t Dim>
    inline H5DataSpace<Dim>& H5DataSpace<Dim>::operator=(const H5DataSpace& obj) {
        Base::operator=(obj);
        return *this;
    }

    template<size_t Dim>
    size_t H5DataSpace<Dim>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return getSimpleExtentNdims();
        else
            return Dim;
    }

    template<size_t Dim>
    size_t H5DataSpace<Dim>::getSize(size_t dim) const noexcept {
        SizeArray sizes(getDim());
        [[maybe_unused]] size_t d = getSimpleExtentDims(sizes.data());
        assert(d == Dim);
        return sizes[dim];
    }

    template<size_t Dim>
    inline H5DataSpace<Dim> H5DataSpace<Dim>::makeDataSpace(const SizeArray& dims) {
        return H5DataSpace(H5::DataSpace(Dim == Utils::Dynamic ? dims.getLength() : Dim, dims.data()));
    }

    template<size_t Dim>
    inline H5DataSpace<Dim> H5DataSpace<Dim>::makeDataSpace(const SizeArray& dims, const SizeArray& maxdims) {
        return H5DataSpace(H5::DataSpace(Dim == Utils::Dynamic ? dims.getLength() : Dim, dims.data(), maxdims.data()));
    }
}
