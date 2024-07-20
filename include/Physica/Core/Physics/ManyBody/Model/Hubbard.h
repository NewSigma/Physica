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

#include "LatticeModel.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim>
    class Hubbard : public LatticeModel<Dim> {
        static_assert(!ScalarType::isComplex, "[Error]: Model param must be real");
        using This = Hubbard<ScalarType, Dim>;
        using Base = LatticeModel<Dim>;
        using HopIndexArray = Utils::Array<Utils::Array<size_t>>;
    public:
        using typename Base::IndexType;
    private:
        ScalarType hoppingT;
        ScalarType repelU;
        HopIndexArray hopIndexArr;
    public:
        Hubbard(Base lattice, ScalarType hoppingT_, ScalarType repelU_);
        Hubbard(const This&) = default;
        Hubbard(This&&) noexcept = default;
        ~Hubbard() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getHoppingT() const noexcept { return hoppingT; }
        [[nodiscard]] ScalarType getRepelU() const noexcept { return repelU; }
        [[nodiscard]] const HopIndexArray& getHopIndexArray() const noexcept { return hopIndexArr; }
    private:
        HopIndexArray makeHopIndexArray();
    };

    template<class ScalarType, unsigned int Dim>
    Hubbard<ScalarType, Dim>::Hubbard(Base lattice, ScalarType hoppingT_, ScalarType repelU_)
            : Base(std::move(lattice)), hoppingT(hoppingT_), repelU(repelU_) {
        if constexpr (Dim > 1)
            hopIndexArr = makeHopIndexArray();
    }

    template<class ScalarType, unsigned int Dim>
    void Hubbard<ScalarType, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        hopIndexArr.swap(obj.hopIndexArr);
    }

    template<class ScalarType, unsigned int Dim>
    typename Hubbard<ScalarType, Dim>::HopIndexArray Hubbard<ScalarType, Dim>::makeHopIndexArray() {
        const auto numSite = Base::getNumSuperCellSite();
        HopIndexArray result(numSite);
        Base::forSiteInLattice([this, numSite, &result](IndexType index) {
            const auto& dims = Base::getDims();
            Utils::Array<size_t> hopTargets{};
            hopTargets.reserve(numSite * Dim * 2);
            for (unsigned int dim = 0; dim < Dim; ++dim) {
                IndexType index1 = index;
                index1[dim] = (index1[dim] + 1) % Base::getSuperSize()[dim];
                hopTargets.append(IndexType::toIndex1D(dims, index1));
            }
            hopTargets.squeeze();
            const auto site = IndexType::toIndex1D(dims, index);
            result[site] = std::move(hopTargets);
        });
        return result;
    }
}
