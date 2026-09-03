/*
 * Copyright 2023-2026 Weibo He.
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

#include "DataSpaceImpl/DataSpaceMixin.h"
#include "DataSpaceImpl/SubDataSpace.h"

namespace Physica {
    template<size_t Dim>
    class H5DataSpace : public DataSpaceMixin<H5DataSpace<Dim>>, public H5ID {
        static_assert(Dim > 0, "[Error]: Invalid dim");
        using Base = DataSpaceMixin<H5DataSpace<Dim>>;
        using This = H5DataSpace<Dim>;
    public:
        using SizeArray = Base::SizeArray;
    private:
        using typename Base::SizeType;

        SizeArray selectedCount{};
        SizeArray selectedStart{};
    public:
        H5DataSpace() = default;
        explicit H5DataSpace(hsize_t size);
        explicit H5DataSpace(SizeArray dims);
        explicit H5DataSpace(H5ID id) noexcept;
        H5DataSpace(SizeArray dims, SizeArray maxdims);
        H5DataSpace(const This&) = default;
        H5DataSpace(This&&) noexcept = default;
        ~H5DataSpace() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start);
        template<size_t SubDim>
        [[nodiscard]] SubDataSpace<This, SubDim> tail(size_t from);
        template<size_t SubDim>
        [[nodiscard]] const SubDataSpace<This, SubDim> tail(size_t from) const;

        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] consteval static size_t getDim() noexcept { return Dim; }
        [[nodiscard]] size_t getSize(size_t dim) const noexcept;
        [[nodiscard]] const auto& getSelectedCount() const noexcept { return selectedCount; }
        [[nodiscard]] const auto& getSelectedStart() const noexcept { return selectedStart; }
        /* Setters */
        void setSize(size_t dim, hsize_t newSize);
        /* Static members */
        [[nodiscard]] constexpr static IdentifierType itype() noexcept { return IdentifierType::Dataspace; }
    };

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(hsize_t size) : H5DataSpace(SizeArray{size}) {
        static_assert(Dim == 1, "[Error]: This constructor applies to 1 dim only");
    }

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(SizeArray dims) : H5ID(H5Screate_simple(Dim, dims.data(), nullptr)) {}

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(H5ID id) noexcept : H5ID(std::move(id)) {
        assert(isa<H5DataSpace>());
    }

    template<size_t Dim>
    H5DataSpace<Dim>::H5DataSpace(SizeArray dims, SizeArray maxdims) : H5ID(H5Screate_simple(Dim, dims.data(), maxdims.data())) {}

    template<size_t Dim>
    H5DataSpace<Dim>& H5DataSpace<Dim>::operator=(This obj) noexcept {
        swap(obj);
        return *this;
    }

    template<size_t Dim>
    void H5DataSpace<Dim>::selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) {
        assert(count.size() == getDim());
        assert(start.size() == getDim());
        selectedCount = count;
        selectedStart = start;
        H5Sselect_hyperslab(getHID(), op, start.data(), nullptr, count.data(), nullptr);
    }

    template<size_t Dim>
    template<size_t SubDim>
    SubDataSpace<H5DataSpace<Dim>, SubDim> H5DataSpace<Dim>::tail(size_t from) {
        return SubDataSpace<H5DataSpace<Dim>, SubDim>(*this, from, getDim());
    }

    template<size_t Dim>
    template<size_t SubDim>
    const SubDataSpace<H5DataSpace<Dim>, SubDim> H5DataSpace<Dim>::tail(size_t from) const {
        return SubDataSpace<H5DataSpace<Dim>, SubDim>(const_cast<This&>(*this), from, getDim());
    }

    template<size_t Dim>
    void H5DataSpace<Dim>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        H5ID::swap(obj);
        selectedCount.swap(obj.selectedCount);
        selectedStart.swap(obj.selectedStart);
    }

    template<size_t Dim>
    size_t H5DataSpace<Dim>::getSize(size_t dim) const noexcept {
        SizeArray sizes{};
        [[maybe_unused]] auto d = H5Sget_simple_extent_dims(getHID(), sizes.data(), nullptr);
        assert(d == getDim());
        return sizes[dim];
    }

    template<size_t Dim>
    void H5DataSpace<Dim>::setSize(size_t dim, hsize_t newSize) {
        SizeArray sizes{};
        [[maybe_unused]] auto d = H5Sget_simple_extent_dims(getHID(), sizes.data(), nullptr);
        assert(d == getDim());
        sizes[dim] = newSize;
        H5Sset_extent_simple(getHID(), Dim, sizes.data(), nullptr);
    }
}

namespace Physica {
    template<size_t S>
    class Traits<H5DataSpace<S>> {
    public:
        constexpr static size_t Dim = S;
    };
}
