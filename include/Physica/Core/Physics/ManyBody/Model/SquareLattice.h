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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/PlainStruct.h"
#include "LatticeModel.h"

namespace Physica {
    template<int Dim, BoundaryCond BC = BoundaryCond::PBC>
    class SquareLattice : public LatticeModel<Dim> {
        using This = SquareLattice<Dim, BC>;
        using Base = LatticeModel<Dim>;
    public:
        using typename Base::DimArray;
        using typename Base::IndexType;
        using ArgVector = std::conditional<BC == BoundaryCond::TBC, DenseVector<float64, Dim>, PlainStruct<void>>::type;
    private:
        struct Hasher {
            std::size_t operator()(const std::pair<int32_t, int32_t>& pair) const noexcept {
                return std::hash<int64_t>{}((static_cast<int64_t>(pair.first) << 32U) | pair.second);
            }
        };

        using NeighborArray = std::conditional<(Dim > 1), Array<Array<size_t>>, PlainStruct<void>>::type;
        using SiteBoundaryMap = std::conditional<BC != BoundaryCond::PBC,
                                                 std::unordered_map<std::pair<int, int>, int, Hasher>,
                                                 PlainStruct<void>>::type;

        [[no_unique_address]] NeighborArray neighbors;
        [[no_unique_address]] NeighborArray nNeighbors;
        [[no_unique_address]] SiteBoundaryMap siteBoundaryMap;
        [[no_unique_address]] ArgVector phaseArgs;
    public:
        SquareLattice() = default;
        SquareLattice(DimArray superSize, size_t numUnitCellSite) requires(BC != BoundaryCond::TBC);
        SquareLattice(DimArray superSize, size_t numUnitCellSite, ArgVector phaseArgs_) requires(BC == BoundaryCond::TBC);
        SquareLattice(const This&) = default;
        SquareLattice(This&&) noexcept = default;
        ~SquareLattice() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void forNeighSites(std::invocable<int, int> auto fn, int site) const;
        void forNNeighSites(std::invocable<int, int> auto fn, int site) const;
        template<Scalar T>
        auto calcPhase() const noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getNeighbors() const noexcept;
        [[nodiscard]] const auto& getNNeighbors() const noexcept;
        [[nodiscard]] const auto& getSiteBoundaryMap() const noexcept;
        [[nodiscard]] const auto& getPhaseArgs() const noexcept;
        /* Setters */
        void setPhaseArgs(ArgVector phaseArgs_) noexcept requires(BC == BoundaryCond::TBC);
    private:
        NeighborArray makeNeighbors() const noexcept;
        NeighborArray makeNNeighbors() const noexcept;
        SiteBoundaryMap makeSiteBoundaryMap() const noexcept;
    };

    template<int Dim, BoundaryCond BC>
    SquareLattice<Dim, BC>::SquareLattice(DimArray superSize, size_t numUnitCellSite) requires(BC != BoundaryCond::TBC)
            : Base(superSize, numUnitCellSite) {
        if constexpr (Dim > 1) {
            neighbors = makeNeighbors();
            nNeighbors = makeNNeighbors();
        }

        if constexpr (BC != BoundaryCond::PBC)
            siteBoundaryMap = makeSiteBoundaryMap();
    }

    template<int Dim, BoundaryCond BC>
    SquareLattice<Dim, BC>::SquareLattice(DimArray superSize, size_t numUnitCellSite, ArgVector phaseArgs_) requires(BC == BoundaryCond::TBC)
            : Base(std::move(superSize), numUnitCellSite) {
        if constexpr (Dim > 1)
            neighbors = makeNeighbors();
        siteBoundaryMap = makeSiteBoundaryMap();
        setPhaseArgs(std::move(phaseArgs_));
    }

    template<int Dim, BoundaryCond BC>
    void SquareLattice<Dim, BC>::forNeighSites(std::invocable<int, int> auto fn, int site) const {
        if constexpr (Dim == 1)
            fn(site, (site + 1) % Base::getNumSuperCellSite());
        else {
            const auto& sites = getNeighbors()[site];
            for (int i : sites)
                fn(site, i);
        }
    }

    template<int Dim, BoundaryCond BC>
    void SquareLattice<Dim, BC>::forNNeighSites(std::invocable<int, int> auto fn, int site) const {
        if constexpr (Dim == 1)
            fn(site, (site + 2) % Base::getNumSuperCellSite());
        else {
            const auto& sites = getNNeighbors()[site];
            for (int i : sites)
                fn(site, i);
        }
    }

    template<int Dim, BoundaryCond BC>
    template<Scalar T>
    auto SquareLattice<Dim, BC>::calcPhase() const noexcept {
        static_assert(!T::isComplex(), "[Error]: Phase arg is real");
        if constexpr (BC == BoundaryCond::TBC) {
            using Tc = T::ComplexType;
            DenseVector<Tc, Dim> result{};
            for (int i = 0; i < Dim; ++i)
                result[i] = Tc::fromPhase(phaseArgs[i]);
            return result;
        }
        else {
            static_assert(BC == BoundaryCond::APBC, "[Error]: Unexpected boundary condition");
            return DenseVector<T, Dim>(Dim, -1);
        }
    }

    template<int Dim, BoundaryCond BC>
    void SquareLattice<Dim, BC>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        neighbors.swap(obj.neighbors);
        siteBoundaryMap.swap(obj.siteBoundaryMap);
        phaseArgs.swap(obj.phaseArgs);
    }

    template<int Dim, BoundaryCond BC>
    const auto& SquareLattice<Dim, BC>::getNeighbors() const noexcept {
        static_assert(Dim > 1, "[Error]: Not available");
        return neighbors;
    }

    template<int Dim, BoundaryCond BC>
    const auto& SquareLattice<Dim, BC>::getNNeighbors() const noexcept {
        static_assert(Dim > 1, "[Error]: Not available");
        return nNeighbors;
    }

    template<int Dim, BoundaryCond BC>
    const auto& SquareLattice<Dim, BC>::getSiteBoundaryMap() const noexcept {
        static_assert(sizeof(siteBoundaryMap) > 1, "[Error]: Not available");
        assert(!siteBoundaryMap.empty());
        return siteBoundaryMap;
    }

    template<int Dim, BoundaryCond BC>
    const auto& SquareLattice<Dim, BC>::getPhaseArgs() const noexcept {
        static_assert(sizeof(phaseArgs) > 1, "[Error]: Not available");
        return phaseArgs;
    }

    template<int Dim, BoundaryCond BC>
    void SquareLattice<Dim, BC>::setPhaseArgs(ArgVector phaseArgs_) noexcept requires(BC == BoundaryCond::TBC) {
        phaseArgs = std::move(phaseArgs_);
    }

    template<int Dim, BoundaryCond BC>
    auto SquareLattice<Dim, BC>::makeNeighbors() const noexcept -> NeighborArray {
        const auto numSite = Base::getNumSuperCellSite();
        NeighborArray result(numSite);
        Base::forSiteInLattice([this, &result](const IndexType index) noexcept {
            const auto dims = Base::calcDims();
            const auto site = IndexType::toIndex1D(dims, index);
            Array<size_t, Dim> hopTargets{};
            for (int dim = 0; dim < Dim; ++dim) {
                IndexType indexN = index;
                indexN[dim] = (indexN[dim] + 1) % Base::getSuperSize()[dim];
                hopTargets[dim] = IndexType::toIndex1D(dims, indexN);
            }
            result[site] = hopTargets;
        });
        return result;
    }

    template<int Dim, BoundaryCond BC>
    auto SquareLattice<Dim, BC>::makeNNeighbors() const noexcept -> NeighborArray {
        static_assert(Dim == 2, "[Error]: Not implemented");
        const auto numSite = Base::getNumSuperCellSite();
        NeighborArray result(numSite);
        Base::forSiteInLattice([this, &result](IndexType index) noexcept {
            const auto dims = Base::calcDims();
            const auto site = IndexType::toIndex1D(dims, index);
            Array<size_t, Dim> hopTargets{};
            index[0] = (index[0] + 1) % dims[0];
            for (int dim = 0; dim < Dim; ++dim) {
                IndexType indexNN = index;
                size_t size = Base::getSuperSize()[1];
                if (dim == 0)
                    indexNN[1] = (index[1] + 1) % size;
                else
                    indexNN[1] = (index[1] + size - 1) % size;
                hopTargets[dim] = IndexType::toIndex1D(dims, indexNN);
            }
            result[site] = hopTargets;
        });
        return result;
    }

    template<int Dim, BoundaryCond BC>
    auto SquareLattice<Dim, BC>::makeSiteBoundaryMap() const noexcept -> SiteBoundaryMap {
        SiteBoundaryMap map{};
        Base::forSiteInLattice([this, &map](const IndexType index) noexcept {
            const auto dims = Base::calcDims();
            const auto site = IndexType::toIndex1D(dims, index);
            for (int dim = 0; dim < Dim; ++dim) {
                bool onBoundary = index[dim] == dims[dim] - 1;
                if (onBoundary) {
                    IndexType index1 = index;
                    index1[dim] = 0;
                    const auto site1 = IndexType::toIndex1D(dims, index1);
                    map[std::make_pair(site, site1)] = dim;
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
