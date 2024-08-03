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

#include "DataSpaceImpl/DataSpaceBase.h"
#include "DataSpaceImpl/SubDataSpace.h"

namespace Physica::Core {
    template<size_t Dim>
    class H5DataSpace : public DataSpaceBase<H5DataSpace<Dim>>, public H5::DataSpace {
        using Base = DataSpaceBase<H5DataSpace<Dim>>;
        using ImplType = H5::DataSpace;
        using This = H5DataSpace<Dim>;
    public:
        using SizeArray = typename Base::SizeArray;
    private:
        using SizeType = typename SizeArray::ValueType;

        SizeArray selectedCount;
        SizeArray selectedStart;
    public:
        explicit H5DataSpace(hsize_t size);
        explicit H5DataSpace(const SizeArray& dims);
        explicit H5DataSpace(std::initializer_list<SizeType> dims);
        H5DataSpace(const SizeArray& dims, const SizeArray& maxdims);
        explicit H5DataSpace(const H5::DataSpace& obj);
        H5DataSpace(const H5DataSpace&) = default;
        H5DataSpace(H5DataSpace&&) noexcept = delete;
        virtual ~H5DataSpace() = default;
        /* Operators */
        inline H5DataSpace& operator=(const H5DataSpace& obj);
        H5DataSpace& operator=(H5DataSpace&&) noexcept = delete;
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start);
        template<size_t SubDim>
        [[nodiscard]] inline SubDataSpace<This, SubDim> tail(size_t from);
        template<size_t SubDim>
        [[nodiscard]] inline const SubDataSpace<This, SubDim> tail(size_t from) const;
        /* Getters */
        [[nodiscard]] const ImplType& asH5Type() const noexcept { return *this; }
        [[nodiscard]] ImplType& asH5Type() noexcept { return *this; }
        [[nodiscard]] inline size_t getDim() const noexcept;
        [[nodiscard]] inline size_t getSize(size_t dim) const noexcept;
        [[nodiscard]] const SizeArray& getSelectedCount() const noexcept { return selectedCount; }
        [[nodiscard]] const SizeArray& getSelectedStart() const noexcept { return selectedStart; }
        [[nodiscard]] bool isValid() const noexcept { return ImplType::isValid(ImplType::getId()); }
        /* Setters */
        void setSize(size_t dim, hsize_t newSize);
    };

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(hsize_t size) : H5DataSpace({size}) {
        static_assert(Dim == 1, "[Error]: This constructor applies to 1 dim only");
    }

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(const SizeArray& dims)
            : H5DataSpace(H5::DataSpace(Dim == Dynamic ? dims.getLength() : Dim, dims.data())) {}

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(std::initializer_list<SizeType> dims) : H5DataSpace(SizeArray(std::move(dims))) {}

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(const SizeArray& dims, const SizeArray& maxdims)
            : H5DataSpace(H5::DataSpace(Dim == Dynamic ? dims.getLength() : Dim, dims.data(), maxdims.data())) {}

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(const H5::DataSpace& obj) : ImplType(obj) {}

    template<size_t Dim>
    void H5DataSpace<Dim>::selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) {
        assert(count.getLength() == getDim());
        assert(start.getLength() == getDim());
        selectedCount = count;
        selectedStart = start;
        ImplType::selectHyperslab(op, count.data(), start.data());
    }

    template<size_t Dim>
    template<size_t SubDim>
    inline SubDataSpace<H5DataSpace<Dim>, SubDim> H5DataSpace<Dim>::tail(size_t from) {
        return SubDataSpace<H5DataSpace<Dim>, SubDim>(*this, from, getDim());
    }

    template<size_t Dim>
    template<size_t SubDim>
    inline const SubDataSpace<H5DataSpace<Dim>, SubDim> H5DataSpace<Dim>::tail(size_t from) const {
        return SubDataSpace<H5DataSpace<Dim>, SubDim>(const_cast<This&>(*this), from, getDim());
    }

    template<size_t Dim>
    inline H5DataSpace<Dim>& H5DataSpace<Dim>::operator=(const H5DataSpace& obj) {
        ImplType::operator=(obj);
        return *this;
    }

    template<size_t Dim>
    inline size_t H5DataSpace<Dim>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return getSimpleExtentNdims();
        else
            return Dim;
    }

    template<size_t Dim>
    inline size_t H5DataSpace<Dim>::getSize(size_t dim) const noexcept {
        SizeArray sizes(getDim());
        [[maybe_unused]] size_t d = getSimpleExtentDims(sizes.data());
        assert(d == getDim());
        return sizes[dim];
    }

    template<size_t Dim>
    void H5DataSpace<Dim>::setSize(size_t dim, hsize_t newSize) {
        SizeArray sizes(getDim());
        [[maybe_unused]] size_t d = getSimpleExtentDims(sizes.data());
        assert(d == getDim());
        sizes[dim] = newSize;
        setExtentSimple(sizes.getLength(), sizes.data());
    }
}

namespace Physica {
    template<size_t S>
    class Traits<Core::H5DataSpace<S>> {
    public:
        constexpr static size_t Dim = S;
    };
}
