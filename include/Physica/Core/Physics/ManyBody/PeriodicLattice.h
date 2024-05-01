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

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    template<unsigned int Dim>
    class PeriodicLattice {
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
    public:
        using DimArray = Utils::Array<unsigned int, Dim>;
        using SiteIndex = Utils::Array<unsigned int, Dim + 1>;
    private:
        DimArray dims;
        unsigned int numSite;
    public:
        PeriodicLattice() = default;
        PeriodicLattice(DimArray dims_, unsigned int numSite_);
        PeriodicLattice(const PeriodicLattice&) = default;
        PeriodicLattice(PeriodicLattice&&) noexcept = default;
        ~PeriodicLattice() = default;
        /* Operators */
        PeriodicLattice& operator=(PeriodicLattice obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor> void forSiteInLattice(Functor func);

        void swap(PeriodicLattice& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] DimArray getDims() const noexcept { return dims; }
        [[nodiscard]] unsigned int getNumUnitCellSite() const noexcept { return numSite; }
        [[nodiscard]] unsigned int getNumSuperCellSite() const noexcept;
    };

    template<unsigned int Dim>
    PeriodicLattice<Dim>::PeriodicLattice(DimArray dims_, unsigned int numSite_) : dims(std::move(dims_)), numSite(numSite_) {}

    template<unsigned int Dim>
    template<class Functor>
    void PeriodicLattice<Dim>::forSiteInLattice(Functor func) {
        if constexpr (Dim == 1) {
            for (unsigned int x = 0; x < dims[0]; ++x)
                for (unsigned int site = 0; site < numSite; ++site)
                    func(SiteIndex{x, site});
        }
        if constexpr (Dim == 2) {
            for (unsigned int x = 0; x < dims[0]; ++x)
                for (unsigned int y = 0; y < dims[1]; ++y)
                    for (unsigned int site = 0; site < numSite; ++site)
                        func(SiteIndex{x, y, site});
        }
        if constexpr (Dim == 3) {
            for (unsigned int x = 0; x < dims[0]; ++x)
                for (unsigned int y = 0; y < dims[1]; ++y)
                    for (unsigned int z = 0; z < dims[2]; ++z)
                        for (unsigned int site = 0; site < numSite; ++site)
                            func(SiteIndex{x, y, z, site});
        }
    }

    template<unsigned int Dim>
    void PeriodicLattice<Dim>::swap(PeriodicLattice& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        dims.swap(obj.dims);
        std::swap(numSite, obj.numSite);
    }

    template<unsigned int Dim>
    unsigned int PeriodicLattice<Dim>::getNumSuperCellSite() const noexcept {
        unsigned int result = getNumUnitCellSite();
        for (auto dim : dims)
            result *= dim;
        return result;
    }
}
