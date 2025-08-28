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

#include <unordered_map>
#include "Physica/PlainStruct.h"
#include "LatticeModel.h"

namespace Physica {
    template<int Dim, BoundaryCond BC = BoundaryCond::PBC>
    class SquareLattice : public LatticeModel<Dim> {
        using This = SquareLattice<Dim, BC>;
        using Base = LatticeModel<Dim>;
    public:
        constexpr static bool UntrivialNearestNeighbor = Dim > 1;
        using typename Base::DimArray;
        using typename Base::IndexType;
    private:
        struct Hash {
            std::size_t operator()(const std::pair<int32_t, int32_t>& pair) const noexcept {
                return std::hash<int64_t>{}((static_cast<int64_t>(pair.first) << 32U) | pair.second);
            }
        };

        using HopIndexArray = std::conditional<UntrivialNearestNeighbor, Array<Array<size_t>>, PlainStruct<void>>::type;
        using SiteBoundaryMap = std::conditional<BC == BoundaryCond::TBC,
                                                 std::unordered_map<std::pair<int, int>, int, Hash>,
                                                 PlainStruct<void>>::type;

        [[no_unique_address]] HopIndexArray hopIndexArr;
        [[no_unique_address]] SiteBoundaryMap siteBoundaryMap;
    public:
        SquareLattice() = default;
        SquareLattice(DimArray superSize, size_t numUnitCellSite);
        SquareLattice(const This&) = default;
        SquareLattice(This&&) noexcept = default;
        ~SquareLattice() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void forNeighSites(std::invocable<int, int> auto fn, int site) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getHopIndexArray() const noexcept { return hopIndexArr; }
        [[nodiscard]] const auto& getSiteBoundaryMap() const noexcept { return siteBoundaryMap; }
    private:
        HopIndexArray makeHopIndexArray() const noexcept;
        SiteBoundaryMap makeSiteBoundaryMap() const noexcept;
    };

    template<int Dim, BoundaryCond BC>
    SquareLattice<Dim, BC>::SquareLattice(DimArray superSize, size_t numUnitCellSite) : Base(superSize, numUnitCellSite) {
        if constexpr (UntrivialNearestNeighbor)
            hopIndexArr = makeHopIndexArray();
        if constexpr (BC == BoundaryCond::TBC)
            siteBoundaryMap = makeSiteBoundaryMap();
    }

    template<int Dim, BoundaryCond BC>
    void SquareLattice<Dim, BC>::forNeighSites(std::invocable<int, int> auto fn, int site) const {
        if constexpr (Dim == 1)
            fn(site, (site + 1) % Base::getNumSuperCellSite());
        else {
            const auto& hopTargets = getHopIndexArray()[site];
            for (int site1 : hopTargets)
                fn(site, site1);
        }
    }

    template<int Dim, BoundaryCond BC>
    void SquareLattice<Dim, BC>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        hopIndexArr.swap(obj.hopIndexArr);
        siteBoundaryMap.swap(obj.siteBoundaryMap);
    }

    template<int Dim, BoundaryCond BC>
    auto SquareLattice<Dim, BC>::makeHopIndexArray() const noexcept -> HopIndexArray {
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

    template<int Dim, BoundaryCond BC>
    auto SquareLattice<Dim, BC>::makeSiteBoundaryMap() const noexcept -> SiteBoundaryMap {
        SiteBoundaryMap map{};
        Base::forSiteInLattice([this, &map](const IndexType index) noexcept {
            const auto& dims = Base::getDims();
            const auto site = IndexType::toIndex1D(dims, index);
            for (int dim = 0; dim < Dim; ++dim) {
                bool onBoundary = index[dim] == dims[dim] - 1;
                if (onBoundary) {
                    IndexType index1 = index;
                    index1[dim] = 0;
                    const auto site1 = IndexType::toIndex1D(dims, index1);
                    map[std::make_pair(site, site1)] = dim;
                    return;
                }
            }
        });
        return map;
    }
}

namespace Physica {
    template<int D, BoundaryCond BC>
    class Traits<SquareLattice<D, BC>> {
    public:
        constexpr static int Dim = D;
        constexpr static int SiteDOF = 4;
    };
}
