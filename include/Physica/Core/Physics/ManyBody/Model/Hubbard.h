/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/PlainStruct.h"
#include "LatticeModel.h"

namespace Physica {
    template<Scalar T, int Dim>
    class Hubbard : public LatticeModel<Dim> {
        static_assert(!T::isComplex, "[Error]: Model param must be real");
        using This = Hubbard<T, Dim>;
        using Base = LatticeModel<Dim>;
    public:
        constexpr static bool UntrivialNearestNeighbor = Dim > 1;
        using typename Base::IndexType;
    private:
        using HopIndexArray = std::conditional<UntrivialNearestNeighbor, Array<Array<size_t>>, PlainStruct<void>>::type;

        T hoppingT;
        T repelU;
        [[no_unique_address]] HopIndexArray hopIndexArr;
    public:
        Hubbard() = default;
        Hubbard(Base lattice, T hoppingT_, T repelU_);
        Hubbard(const This&) = default;
        Hubbard(This&&) noexcept = default;
        ~Hubbard() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getHoppingT() const noexcept { return hoppingT; }
        [[nodiscard]] T getRepelU() const noexcept { return repelU; }
        [[nodiscard]] const auto& getHopIndexArray() const noexcept { return hopIndexArr; }
    private:
        HopIndexArray makeHopIndexArray() noexcept;
    };

    template<Scalar T, int Dim>
    Hubbard<T, Dim>::Hubbard(Base lattice, T hoppingT_, T repelU_)
            : Base(std::move(lattice)), hoppingT(hoppingT_), repelU(repelU_) {
        if constexpr (UntrivialNearestNeighbor)
            hopIndexArr = makeHopIndexArray();
    }

    template<Scalar T, int Dim>
    void Hubbard<T, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        hopIndexArr.swap(obj.hopIndexArr);
    }

    template<Scalar T, int Dim>
    auto Hubbard<T, Dim>::makeHopIndexArray() noexcept -> HopIndexArray {
        const auto numSite = Base::getNumSuperCellSite();
        HopIndexArray result(numSite);
        Base::forSiteInLattice([this, numSite, &result](const IndexType index) noexcept {
            const auto& dims = Base::getDims();
            Array<size_t> hopTargets{};
            hopTargets.reserve(numSite * Dim * 2);
            for (int dim = 0; dim < Dim; ++dim) {
                const bool isMultiCount = dims[dim] == 2 && index[dim] == 1;
                if (isMultiCount)
                    continue;

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

namespace Physica {
    template<Scalar T, int Dim_>
    class Traits<Hubbard<T, Dim_>> {
    public:
        constexpr static int Dim = Dim_;
        constexpr static int SiteDOF = 4;
    };
}
