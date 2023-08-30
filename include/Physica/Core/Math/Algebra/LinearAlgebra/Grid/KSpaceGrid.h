/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Physics/SolidState/PeriodicCell.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "RSpaceGrid.h"

namespace Physica::Core {
    template<class GridType> class KSpaceGrid;

    namespace Internal {
        template<class GridType>
        class Traits<KSpaceGrid<GridType>> {

        public:
            using GridScalar = typename GridType::ScalarType;
            using ScalarType = typename GridScalar::ComplexType;
            constexpr static bool isComplex = GridScalar::isComplex;
        };
    }

    template<class GridType>
    class KSpaceGrid : public LValueGrid<KSpaceGrid<GridType>> {
        using This = KSpaceGrid<GridType>;
        using Base = LValueGrid<This>;
        using Traits = Internal::Traits<This>;
        using GridScalar = typename Traits::GridScalar;
        constexpr static bool isComplex = Traits::isComplex;
    public:
        using ScalarType = typename Base::ScalarType;
        using RealType = typename ScalarType::RealType;
        using Index3D = typename Base::Index3D;
        using VectorType = typename Base::VectorType;
        using LatticeMatrix = typename Base::LatticeMatrix;
    private:
        RSpaceGrid<ScalarType> values;
        size_t rSpaceSizeZ;
    public:
        KSpaceGrid() = default;
        template<class... Args>
        KSpaceGrid(Index3D rSpaceDim, Args... args);
        KSpaceGrid(RSpaceGrid<ScalarType> values_, size_t rSpaceSizeZ_);
        KSpaceGrid(const KSpaceGrid&) = default;
        KSpaceGrid(KSpaceGrid&&) noexcept = default;
        ~KSpaceGrid() = default;
        /* Operators */
        using Base::operator();
        KSpaceGrid& operator=(KSpaceGrid grid) noexcept;
        [[nodiscard]] inline ScalarType& operator()(size_t x, size_t y, size_t z);
        [[nodiscard]] inline const ScalarType& operator()(size_t x, size_t y, size_t z) const;
        template<class U>
        friend std::ostream& operator<<(std::ostream& os, const KSpaceGrid<U>& grid);
        template<class U>
        friend std::istream& operator>>(std::istream& is, KSpaceGrid<U>& grid);
        /* Operations */
        inline void swap(KSpaceGrid& grid) noexcept;
        /* Getters */
        [[nodiscard]] RSpaceGrid<ScalarType>& getValues() { return values; }
        [[nodiscard]] const RSpaceGrid<ScalarType>& getValues() const { return values; }
        [[nodiscard]] size_t getDimX() const noexcept { return values.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return values.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return rSpaceSizeZ; }
        using Base::flatten;
        /* Static members */
        static typename Base::Index3D makeGridDim(RealType cutEnergy, const LatticeMatrix& repCell);
        static KSpaceGrid<GridType> makeGrid(RealType cutEnergy, const LatticeMatrix& repCell);
        template<class Functor>
        inline static void forKInGrid(const KSpaceGrid& grid, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        inline static void forReducedKInGrid(const KSpaceGrid& grid, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        inline static void forKIndexInGrid(const KSpaceGrid& grid, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        inline static void forReducedKIndexInGrid(const KSpaceGrid& grid, const LatticeMatrix& repLatt, Functor func);
    private:
        using Base::forPointInGrid;
        using Base::forPointIndexInGrid;
        using Base::forIndexInGrid;
    };

    template<class GridType>
    template<class... Args>
    KSpaceGrid<GridType>::KSpaceGrid(Index3D rSpaceDim, Args... args)
            : values(FFT<GridScalar, 3>::rSizeToKSize(rSpaceDim), std::forward<Args>(args)...)
            , rSpaceSizeZ(rSpaceDim[2]) {}

    template<class GridType>
    KSpaceGrid<GridType>::KSpaceGrid(RSpaceGrid<ScalarType> values_, size_t rSpaceSizeZ_)
            : values(std::move(values_)), rSpaceSizeZ(rSpaceSizeZ_) {
        assert((FFT<GridScalar, 1>::rSizeToKSize(rSpaceSizeZ) == values.getDimZ()) && "[Error]: Incompatible size");
    }

    template<class GridType>
    KSpaceGrid<GridType>& KSpaceGrid<GridType>::operator=(KSpaceGrid<GridType> grid) noexcept {
        swap(grid);
        return *this;
    }

    template<class GridType>
    inline typename KSpaceGrid<GridType>::ScalarType& KSpaceGrid<GridType>::operator()(size_t x, size_t y, size_t z) {
        assert(x < getDimX());
        assert(y < getDimY());
        assert(z < getDimZ());
        if constexpr (isComplex)
            return values(x, y, z);
        else {
            const bool isOverflow = z >= values.getDimZ();
            const size_t x1 = (isOverflow && x != 0) ? (getDimX() - x) : x;
            const size_t y1 = (isOverflow && y != 0) ? (getDimY() - y) : y;
            const size_t z1 = getDimZ() - z;
            return values(x1, y1, z1);
        }
    }

    template<class GridType>
    inline const typename KSpaceGrid<GridType>::ScalarType& KSpaceGrid<GridType>::operator()(size_t x, size_t y, size_t z) const {
        return const_cast<This&>(*this).operator()(x, y, z);
    }

    template<class GridType>
    std::ostream& operator<<(std::ostream& os, const KSpaceGrid<GridType>& grid) {
        os << grid.values;
        return os;
    }
    template<class GridType>
    std::istream& operator>>(std::istream& is, KSpaceGrid<GridType>& grid) {
        is >> grid.values;
        return is;
    }

    template<class GridType>
    void KSpaceGrid<GridType>::swap(KSpaceGrid& grid) noexcept {
        values.swap(grid.values);
        std::swap(rSpaceSizeZ, grid.rSpaceSizeZ);
    }

    template<class GridType>
    typename KSpaceGrid<GridType>::Index3D KSpaceGrid<GridType>::makeGridDim(RealType cutEnergy, const LatticeMatrix& repCell) {
        constexpr double factor = 2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck;
        const RealType maxWaveVec = sqrt(RealType(factor) * cutEnergy);
        const auto range = PeriodicCell<RealType, 3>::estimateRange(repCell, maxWaveVec);
        return {size_t(2 * range[0] + 1), size_t(2 * range[1] + 1), size_t(range[2] + 1)};
    }

    template<class GridType>
    KSpaceGrid<GridType> KSpaceGrid<GridType>::makeGrid(RealType cutEnergy, const LatticeMatrix& repCell) {
        return KSpaceGrid<GridType>(makeGridDim(cutEnergy, repCell));
    }

    template<class GridType>
    template<class Functor>
    inline void KSpaceGrid<GridType>::forKInGrid(const KSpaceGrid<GridType>& grid, const LatticeMatrix& repLatt, Functor func) {
        RSpaceGrid<ScalarType>::template forPointInGrid<false, Functor>(grid.getDim(), repLatt, func);
    }

    template<class GridType>
    template<class Functor>
    inline void KSpaceGrid<GridType>::forReducedKInGrid(const KSpaceGrid<GridType>& grid, const LatticeMatrix& repLatt, Functor func) {
        LatticeMatrix copy = repLatt;
        auto row = copy.row(2);
        row *= RealType(grid.values.getDimZ()) / RealType(grid.getDimZ());
        RSpaceGrid<ScalarType>::template forPointInGrid<false, Functor>(grid.values.getDim(), copy, func);
    }

    template<class GridType>
    template<class Functor>
    inline void KSpaceGrid<GridType>::forKIndexInGrid(const KSpaceGrid<GridType>& grid, const LatticeMatrix& repLatt, Functor func) {
        RSpaceGrid<ScalarType>::template forPointIndexInGrid<false, Functor>(grid.getDim(), repLatt, func);
    }

    template<class GridType>
    template<class Functor>
    inline void KSpaceGrid<GridType>::forReducedKIndexInGrid(const KSpaceGrid<GridType>& grid, const LatticeMatrix& repLatt, Functor func) {
        LatticeMatrix copy = repLatt;
        auto row = copy.row(2);
        row *= RealType(grid.values.getDimZ()) / RealType(grid.getDimZ());
        RSpaceGrid<ScalarType>::template forPointIndexInGrid<false, Functor>(grid.values.getDim(), repLatt, func);
    }
}
