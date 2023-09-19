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

#include "SpinPair.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::Core {
    template<class ScalarType, bool IsSpinPolarized>
    class DensityGrid {
        using VectorType = typename RSpaceGrid<ScalarType>::VectorType;
        using LatticeMatrix = typename RSpaceGrid<ScalarType>::LatticeMatrix;
        using GridType = RSpaceGrid<ScalarType>;
    public:
        using Index3D = typename RSpaceGrid<ScalarType>::Index3D;
        using BasisType = PlainWaveBasis<ScalarType>;
        using KSOrbitArray = Utils::Array<SpinPair<BasisType, IsSpinPolarized>>;
    private:
        SpinPair<GridType, IsSpinPolarized> densityPair;
    public:
        DensityGrid() = default;
        DensityGrid(Index3D dim);
        DensityGrid(const DensityGrid&) = default;
        DensityGrid(DensityGrid&&) noexcept = default;
        ~DensityGrid() = default;
        /* Oprators */
        DensityGrid& operator=(DensityGrid obj) noexcept;
        [[nodiscard]] ScalarType operator()(SpinState spin, VectorType pos) const;
        /* Operations */
        void initDensity(ScalarType averageRho);
        inline void initDensity(const DensityGrid& rho);

        void resize(size_t x, size_t y, size_t z) { resize({x, y, z}); }
        void resize(Index3D index);
        void fit(const DensityGrid& rho);
        void swap(DensityGrid& obj) noexcept;
        /* Getters */
        [[nodiscard]] const GridType& getTotalDensity() const noexcept { return densityPair[SpinState::Up]; }
        [[nodiscard]] GridType& getTotalDensity() noexcept { return densityPair[SpinState::Up]; }
        [[nodiscard]] const GridType& getPolarDensity() const noexcept { return densityPair[SpinState::Down]; }
        [[nodiscard]] GridType& getPolarDensity() noexcept { return densityPair[SpinState::Down]; }
    };

    template<class ScalarType, bool IsSpinPolarized>
    DensityGrid<ScalarType, IsSpinPolarized>::DensityGrid(Index3D dim) : densityPair(dim) {}

    template<class ScalarType, bool IsSpinPolarized>
    DensityGrid<ScalarType, IsSpinPolarized>& DensityGrid<ScalarType, IsSpinPolarized>::operator=(DensityGrid obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, bool IsSpinPolarized>
    ScalarType DensityGrid<ScalarType, IsSpinPolarized>::operator()(SpinState spin, VectorType pos) const {
        using PosScalarType = typename VectorType::ScalarType;
        assert(PosScalarType(0) <= pos[0] && pos[0] <= PosScalarType(1));
        assert(PosScalarType(0) <= pos[1] && pos[1] <= PosScalarType(1));
        assert(PosScalarType(0) <= pos[2] && pos[2] <= PosScalarType(1));
        const size_t dimX = getTotalDensity().getDimX();
        const size_t dimY = getTotalDensity().getDimY();
        const size_t dimZ = getTotalDensity().getDimZ();
        const size_t nx1 = double(PosScalarType(dimX) * pos[0]);
        const size_t ny1 = double(PosScalarType(dimY) * pos[1]);
        const size_t nz1 = double(PosScalarType(dimZ) * pos[2]);
        const size_t nx2 = (nx1 + 1) % dimX;
        const size_t ny2 = (ny1 + 1) % dimY;
        const size_t nz2 = (nz1 + 1) % dimZ;

        const ScalarType deltaX = reciprocal(ScalarType(dimX));
        const ScalarType deltaY = reciprocal(ScalarType(dimY));
        const ScalarType deltaZ = reciprocal(ScalarType(dimZ));
        const ScalarType x1 = deltaX * ScalarType(nx1);
        const ScalarType y1 = deltaY * ScalarType(ny1);
        const ScalarType z1 = deltaZ * ScalarType(nz1);
        const ScalarType factorX2 = (pos[0] - x1) / deltaX;
        const ScalarType factorY2 = (pos[1] - y1) / deltaY;
        const ScalarType factorZ2 = (pos[2] - z1) / deltaZ;
        const ScalarType factorX1 = ScalarType(1) - factorX1;
        const ScalarType factorY1 = ScalarType(1) - factorY1;
        const ScalarType factorZ1 = ScalarType(1) - factorZ1;
        const auto& grid = densityPair[spin];
        return grid(nx1, ny1, nz1) * (factorX1 * factorY1 * factorZ1)
             + grid(nx2, ny1, nz1) * (factorX2 * factorY1 * factorZ1)
             + grid(nx1, ny2, nz1) * (factorX1 * factorY2 * factorZ1)
             + grid(nx1, ny1, nz2) * (factorX1 * factorY1 * factorZ2)
             + grid(nx2, ny2, nz1) * (factorX2 * factorY2 * factorZ1)
             + grid(nx2, ny1, nz2) * (factorX2 * factorY1 * factorZ2)
             + grid(nx1, ny2, nz2) * (factorX1 * factorY2 * factorZ2)
             + grid(nx2, ny2, nz2) * (factorX2 * factorY2 * factorZ2);
    }

    template<class ScalarType, bool IsSpinPolarized>
    void DensityGrid<ScalarType, IsSpinPolarized>::resize(Index3D dim) {
        getTotalDensity().resize(dim);
        if constexpr (IsSpinPolarized)
            getPolarDensity().resize(dim);
    }

    template<class ScalarType, bool IsSpinPolarized>
    void DensityGrid<ScalarType, IsSpinPolarized>::initDensity(ScalarType averageRho) {
        {
            auto rho = getTotalDensity().flatten();
            rho = averageRho;
        }
        if constexpr (IsSpinPolarized) {
            auto zeta = getPolarDensity().flatten();
            zeta = ScalarType(0);
        }
    }

    template<class ScalarType, bool IsSpinPolarized>
    inline void DensityGrid<ScalarType, IsSpinPolarized>::initDensity(const DensityGrid& rho) {
        fit(rho);
    }

    template<class ScalarType, bool IsSpinPolarized>
    void DensityGrid<ScalarType, IsSpinPolarized>::fit(const DensityGrid& rho) {
        const LatticeMatrix latt = LatticeMatrix::unitMatrix(3);
        auto kernel = [this, &rho](VectorType pos, Index3D index) {
            {
                auto& rho_up = getTotalDensity();
                rho_up(index) = rho(SpinState::Up, pos);
            }
            if constexpr (IsSpinPolarized) {
                auto& rho_down = getPolarDensity();
                rho_down(index) = rho(SpinState::Down, pos);
            }
        };
        RSpaceGrid<ScalarType>::template forPointIndexInGrid<true, decltype(kernel)>(getTotalDensity(), latt, kernel);
        {
            auto rho_up_new = getTotalDensity().flatten();
            const auto& rho_up_old = rho[SpinState::Up].flatten();
            const ScalarType factor = mean(rho_up_old) / mean(rho_up_new);
            rho_up_new *= factor;
        }
        if constexpr (IsSpinPolarized) {
            auto rho_down_new = getPolarDensity().flatten();
            const auto& rho_down_old = rho[SpinState::Down].flatten();
            const ScalarType factor = mean(rho_down_old) / mean(rho_down_new);
            rho_down_new *= factor;
        }
    }

    template<class ScalarType, bool IsSpinPolarized>
    void DensityGrid<ScalarType, IsSpinPolarized>::swap(DensityGrid& obj) noexcept {
        densityPair.swap(obj.densityPair);
    }

    template<class ScalarType, bool IsSpinPolarized>
    std::ostream& operator<<(std::ostream& os, const DensityGrid<ScalarType, IsSpinPolarized>& grid) {
        os << grid.getTotalDensity();
        if (IsSpinPolarized)
            os << grid.getPolarDensity();
        return os;
    }

    template<class ScalarType, bool IsSpinPolarized>
    std::istream& operator>>(std::istream& is, DensityGrid<ScalarType, IsSpinPolarized>& grid) {
        is >> grid.getTotalDensity();
        if (IsSpinPolarized)
            is >> grid.getPolarDensity();
        return is;
    }
}
