/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h>
#include <Physica/Core/Math/Statistics/NumCharacter.h>
#include "Basis/PlainWaveBasis.h"
#include "SpinPair.h"

namespace Physica::Core {
    template<Scalar T, bool IsSpinPolarized>
    class DensityGrid {
        using LatticeMatrix = typename PeriodicCell<T, 3>::LatticeMatrix;
        using GridType = RSpaceGrid<T>;
    public:
        using Index3D = typename RSpaceGrid<T>::Index3D;
        using BasisType = PlainWaveBasis<T>;
        using KSOrbitArray = Array<SpinPair<BasisType, IsSpinPolarized>>;
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
        [[nodiscard]] T operator()(SpinState spin, Vector3D<T> pos) const;
        /* Operations */
        void initDensity(T averageRho);
        inline void initDensity(const DensityGrid& rho);

        void resize(size_t x, size_t y, size_t z) { resize({x, y, z}); }
        void resize(Index3D index);
        void fit(const DensityGrid& rho);
        void swap(DensityGrid& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const GridType& getTotalDensity() const noexcept { return densityPair[SpinState::Up]; }
        [[nodiscard]] GridType& getTotalDensity() noexcept { return densityPair[SpinState::Up]; }
        [[nodiscard]] const GridType& getPolarDensity() const noexcept { return densityPair[SpinState::Down]; }
        [[nodiscard]] GridType& getPolarDensity() noexcept { return densityPair[SpinState::Down]; }
    };

    template<Scalar T, bool IsSpinPolarized>
    DensityGrid<T, IsSpinPolarized>::DensityGrid(Index3D dim) : densityPair(dim) {}

    template<Scalar T, bool IsSpinPolarized>
    DensityGrid<T, IsSpinPolarized>& DensityGrid<T, IsSpinPolarized>::operator=(DensityGrid obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T, bool IsSpinPolarized>
    T DensityGrid<T, IsSpinPolarized>::operator()(SpinState spin, Vector3D<T> pos) const {
        assert(T(0) <= pos[0] && pos[0] <= T(1));
        assert(T(0) <= pos[1] && pos[1] <= T(1));
        assert(T(0) <= pos[2] && pos[2] <= T(1));
        const size_t dimX = getTotalDensity().getDimX();
        const size_t dimY = getTotalDensity().getDimY();
        const size_t dimZ = getTotalDensity().getDimZ();
        const size_t nx1 = double(T(dimX) * pos[0]);
        const size_t ny1 = double(T(dimY) * pos[1]);
        const size_t nz1 = double(T(dimZ) * pos[2]);
        const size_t nx2 = (nx1 + 1) % dimX;
        const size_t ny2 = (ny1 + 1) % dimY;
        const size_t nz2 = (nz1 + 1) % dimZ;

        const T deltaX = reciprocal(T(dimX));
        const T deltaY = reciprocal(T(dimY));
        const T deltaZ = reciprocal(T(dimZ));
        const T x1 = deltaX * T(nx1);
        const T y1 = deltaY * T(ny1);
        const T z1 = deltaZ * T(nz1);
        const T factorX2 = (pos[0] - x1) / deltaX;
        const T factorY2 = (pos[1] - y1) / deltaY;
        const T factorZ2 = (pos[2] - z1) / deltaZ;
        const T factorX1 = T(1) - factorX2;
        const T factorY1 = T(1) - factorY2;
        const T factorZ1 = T(1) - factorZ2;
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

    template<Scalar T, bool IsSpinPolarized>
    void DensityGrid<T, IsSpinPolarized>::resize(Index3D dim) {
        getTotalDensity().resize(dim);
        if constexpr (IsSpinPolarized)
            getPolarDensity().resize(dim);
    }

    template<Scalar T, bool IsSpinPolarized>
    void DensityGrid<T, IsSpinPolarized>::initDensity(T averageRho) {
        {
            auto rho = getTotalDensity().flatten();
            rho = averageRho;
        }
        if constexpr (IsSpinPolarized) {
            auto zeta = getPolarDensity().flatten();
            zeta = T(0);
        }
    }

    template<Scalar T, bool IsSpinPolarized>
    inline void DensityGrid<T, IsSpinPolarized>::initDensity(const DensityGrid& rho) {
        fit(rho);
    }

    template<Scalar T, bool IsSpinPolarized>
    void DensityGrid<T, IsSpinPolarized>::fit(const DensityGrid& rho) {
        const LatticeMatrix latt = LatticeMatrix::unitMatrix(3);
        auto kernel = [this, &rho](Vector3D<T> pos, Index3D index) {
            {
                auto& rho_up = getTotalDensity();
                rho_up(index) = rho(SpinState::Up, pos);
            }
            if constexpr (IsSpinPolarized) {
                auto& rho_down = getPolarDensity();
                rho_down(index) = rho(SpinState::Down, pos);
            }
        };
        RSpaceGrid<T>::template forPointIndexInGrid<T, true, decltype(kernel)>(getTotalDensity(), latt, kernel);
        {
            auto rho_up_new = getTotalDensity().flatten();
            const auto& rho_up_old = rho.getTotalDensity().flatten();
            const T factor = mean(rho_up_old) / mean(rho_up_new);
            rho_up_new *= factor;
        }
        if constexpr (IsSpinPolarized) {
            auto rho_down_new = getPolarDensity().flatten();
            const auto& rho_down_old = rho.getPolarDensity().flatten();
            const T factor = mean(rho_down_old) / mean(rho_down_new);
            rho_down_new *= factor;
        }
    }

    template<Scalar T, bool IsSpinPolarized>
    void DensityGrid<T, IsSpinPolarized>::swap(DensityGrid& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        densityPair.swap(obj.densityPair);
    }

    template<Scalar T, bool IsSpinPolarized>
    std::ostream& operator<<(std::ostream& os, const DensityGrid<T, IsSpinPolarized>& grid) {
        os << grid.getTotalDensity();
        if (IsSpinPolarized)
            os << grid.getPolarDensity();
        return os;
    }

    template<Scalar T, bool IsSpinPolarized>
    std::istream& operator>>(std::istream& is, DensityGrid<T, IsSpinPolarized>& grid) {
        is >> grid.getTotalDensity();
        if (IsSpinPolarized)
            is >> grid.getPolarDensity();
        return is;
    }
}
