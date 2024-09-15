/*
 * Copyright 2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>
#include "SiteIndex.h"

namespace Physica::Core {
    template<unsigned int Dim>
    class LatticeModel {
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        using This = LatticeModel<Dim>;
    public:
        using DimArray = Array<size_t, Dim>;
        using IndexType = SiteIndex<Dim>;
    private:
        DimArray superSize;
        size_t numUnitCellSite;
    public:
        LatticeModel() = default;
        LatticeModel(DimArray superSize_, size_t numUnitCellSite_);
        LatticeModel(const This&) = default;
        LatticeModel(This&&) noexcept = default;
        ~LatticeModel() = default;
        /* Operators */
        LatticeModel& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor> void forSiteInLattice(Functor func) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const DimArray& getSuperSize() const noexcept { return superSize; }
        [[nodiscard]] size_t getNumCell() const noexcept;
        [[nodiscard]] size_t getNumUnitCellSite() const noexcept { return numUnitCellSite; }
        [[nodiscard]] size_t getNumSuperCellSite() const noexcept { return getNumUnitCellSite() * getNumCell(); }
        [[nodiscard]] IndexType getDims() const noexcept;
    };

    template<unsigned int Dim>
    LatticeModel<Dim>::LatticeModel(DimArray superSize_, size_t numUnitCellSite_)
            : superSize(std::move(superSize_)), numUnitCellSite(numUnitCellSite_) {}

    template<unsigned int Dim>
    template<class Functor>
    void LatticeModel<Dim>::forSiteInLattice(Functor func) const {
        if constexpr (Dim == 1) {
            for (size_t x = 0; x < superSize[0]; ++x)
                for (size_t site = 0; site < numUnitCellSite; ++site)
                    func(IndexType{x, site});
        }
        if constexpr (Dim == 2) {
            for (size_t x = 0; x < superSize[0]; ++x)
                for (size_t y = 0; y < superSize[1]; ++y)
                    for (size_t site = 0; site < numUnitCellSite; ++site)
                        func(IndexType{x, y, site});
        }
        if constexpr (Dim == 3) {
            for (size_t x = 0; x < superSize[0]; ++x)
                for (size_t y = 0; y < superSize[1]; ++y)
                    for (size_t z = 0; z < superSize[2]; ++z)
                        for (size_t site = 0; site < numUnitCellSite; ++site)
                            func(IndexType{x, y, z, site});
        }
    }

    template<unsigned int Dim>
    void LatticeModel<Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        superSize.swap(obj.superSize);
        std::swap(numUnitCellSite, obj.numUnitCellSite);
    }

    template<unsigned int Dim>
    size_t LatticeModel<Dim>::getNumCell() const noexcept {
        size_t result = 1;
        for (auto size : superSize)
            result *= size;
        return result;
    }

    template<unsigned int Dim>
    typename LatticeModel<Dim>::IndexType LatticeModel<Dim>::getDims() const noexcept {
        if constexpr (Dim == 1)
            return IndexType{superSize[0], numUnitCellSite};
        if constexpr (Dim == 2)
            return IndexType{superSize[0], superSize[1], numUnitCellSite};
        if constexpr (Dim == 3)
            return IndexType{superSize[0], superSize[1], superSize[2], numUnitCellSite};
    }
}
