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
#include "Physica/Core/Physics/Container/RSpaceGrid.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::Core {
    template<class ScalarType, bool isSpinPolarized>
    class DensityGrid : public SpinPair<RSpaceGrid<ScalarType>, isSpinPolarized> {
        using Base = SpinPair<RSpaceGrid<ScalarType>, isSpinPolarized>;
        using VectorType = typename RSpaceGrid<ScalarType>::VectorType;
        using LatticeMatrix = typename RSpaceGrid<ScalarType>::LatticeMatrix;
        using Index3D = typename RSpaceGrid<ScalarType>::Index3D;
    public:
        using Base::Base;
        DensityGrid(const DensityGrid&) = default;
        DensityGrid(DensityGrid&&) noexcept = default;
        ~DensityGrid() = default;
        /* Oprators */
        DensityGrid& operator=(DensityGrid rho) noexcept;
        [[nodiscard]] ScalarType operator()(SpinState spin, VectorType pos) const;
        /* Operations */
        void resize(size_t x, size_t y, size_t z);
        void fit(const DensityGrid& rho);
        void swap(DensityGrid& rho) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDimX() const noexcept { return Base::operator[](SpinState::Up).getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return Base::operator[](SpinState::Up).getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return Base::operator[](SpinState::Up).getDimZ(); }
    };

    template<class ScalarType, bool isSpinPolarized>
    DensityGrid<ScalarType, isSpinPolarized>& DensityGrid<ScalarType, isSpinPolarized>::operator=(DensityGrid rho) noexcept {
        swap(rho);
        return *this;
    }

    template<class ScalarType, bool isSpinPolarized>
    void DensityGrid<ScalarType, isSpinPolarized>::resize(size_t x, size_t y, size_t z) {
        Base::operator[](SpinState::Up).resize(x, y, z);
        if constexpr (isSpinPolarized)
            Base::operator[](SpinState::Down).resize(x, y, z);;
    }

    template<class ScalarType, bool isSpinPolarized>
    ScalarType DensityGrid<ScalarType, isSpinPolarized>::operator()(SpinState spin, VectorType pos) const {
        assert(ScalarType(0) <= pos[0] && pos[0] <= ScalarType(1));
        assert(ScalarType(0) <= pos[1] && pos[1] <= ScalarType(1));
        assert(ScalarType(0) <= pos[2] && pos[2] <= ScalarType(1));
        const size_t nx1 = double(ScalarType(getDimX()) * pos[0]);
        const size_t ny1 = double(ScalarType(getDimY()) * pos[1]);
        const size_t nz1 = double(ScalarType(getDimZ()) * pos[2]);
        const size_t nx2 = (nx1 + 1) % getDimX();
        const size_t ny2 = (ny1 + 1) % getDimY();
        const size_t nz2 = (nz1 + 1) % getDimZ();

        const ScalarType deltaX = reciprocal(ScalarType(getDimX()));
        const ScalarType deltaY = reciprocal(ScalarType(getDimY()));
        const ScalarType deltaZ = reciprocal(ScalarType(getDimZ()));
        const ScalarType x1 = deltaX * ScalarType(nx1);
        const ScalarType y1 = deltaY * ScalarType(ny1);
        const ScalarType z1 = deltaZ * ScalarType(nz1);
        const ScalarType factorX2 = (pos[0] - x1) / deltaX;
        const ScalarType factorY2 = (pos[1] - y1) / deltaY;
        const ScalarType factorZ2 = (pos[2] - z1) / deltaZ;
        const ScalarType factorX1 = ScalarType(1) - factorX1;
        const ScalarType factorY1 = ScalarType(1) - factorY1;
        const ScalarType factorZ1 = ScalarType(1) - factorZ1;
        const auto& grid = Base::operator[](spin);
        return grid(nx1, ny1, nz1) * (factorX1 * factorY1 * factorZ1)
             + grid(nx2, ny1, nz1) * (factorX2 * factorY1 * factorZ1)
             + grid(nx1, ny2, nz1) * (factorX1 * factorY2 * factorZ1)
             + grid(nx1, ny1, nz2) * (factorX1 * factorY1 * factorZ2)
             + grid(nx2, ny2, nz1) * (factorX2 * factorY2 * factorZ1)
             + grid(nx2, ny1, nz2) * (factorX2 * factorY1 * factorZ2)
             + grid(nx1, ny2, nz2) * (factorX1 * factorY2 * factorZ2)
             + grid(nx2, ny2, nz2) * (factorX2 * factorY2 * factorZ2);
    }

    template<class ScalarType, bool isSpinPolarized>
    void DensityGrid<ScalarType, isSpinPolarized>::fit(const DensityGrid& rho) {
        const LatticeMatrix latt = LatticeMatrix::unitMatrix(3);
        RSpaceGrid<ScalarType>::forPointIndexInGrid(Base::operator[](SpinState::Up), latt,
            [this, &rho](VectorType pos, Index3D index) {
                {
                    auto& rho_up = Base::operator[](SpinState::Up);
                    rho_up(index) = rho(SpinState::Up, pos);
                }
                if constexpr (isSpinPolarized) {
                    auto& rho_down = Base::operator[](SpinState::Down);
                    rho_down(index) = rho(SpinState::Down, pos);
                }
            });
        {
            auto& rho_up_new = Base::operator[](SpinState::Up).asVector();
            const auto& rho_up_old = rho[SpinState::Up].asVector();
            const ScalarType factor = mean(rho_up_old) / mean(rho_up_new);
            rho_up_new *= factor;
        }
        if constexpr (isSpinPolarized) {
            auto& rho_down_new = Base::operator[](SpinState::Down).asVector();
            const auto& rho_down_old = rho[SpinState::Down].asVector();
            const ScalarType factor = mean(rho_down_old) / mean(rho_down_new);
            rho_down_new *= factor;
        }
    }

    template<class ScalarType, bool isSpinPolarized>
    void DensityGrid<ScalarType, isSpinPolarized>::swap(DensityGrid& rho) noexcept {
        Base::swap(rho);
    }

    template<class ScalarType, bool isSpinPolarized>
    std::ostream& operator<<(std::ostream& os, const DensityGrid<ScalarType, isSpinPolarized>& grid) {
        os << grid[SpinState::Up];
        if (isSpinPolarized)
            os << grid[SpinState::Down];
        return os;
    }

    template<class ScalarType, bool isSpinPolarized>
    std::istream& operator>>(std::istream& is, DensityGrid<ScalarType, isSpinPolarized>& grid) {
        is >> grid[SpinState::Up];
        if (isSpinPolarized)
            is >> grid[SpinState::Down];
        return is;
    }
}
