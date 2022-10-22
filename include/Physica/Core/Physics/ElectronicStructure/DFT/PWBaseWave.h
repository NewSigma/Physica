/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "KSpaceGrid.h"

namespace Physica::Core {
    /**
     * \class PWBaseWave
     * Plain wave base wave function
     */
    template<class ScalarType>
    class PWBaseWave {
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using ComplexType = ComplexScalar<ScalarType>;
        using GridType = KSpaceGrid<ComplexType>;

        GridType grid;
        LatticeMatrix repCell;
    public:
        PWBaseWave(ScalarType cutEnergy, LatticeMatrix repCell_);
        /* Operators */
        [[nodiscard]] ComplexType operator()(Vector<ScalarType, 3> r) const;
        template<class VectorType>
        PWBaseWave& operator=(const RValueVector<VectorType>& newCoeffs);
        /* Getters */
        [[nodiscard]] size_t getNumPlaneWave() const noexcept { return grid.getSize(); }
        [[nodiscard]] ssize_t getDimX() const noexcept { return grid.getDimX(); }
        [[nodiscard]] ssize_t getDimY() const noexcept { return grid.getDimY(); }
        [[nodiscard]] ssize_t getDimZ() const noexcept { return grid.getDimZ(); }
        [[nodiscard]] typename GridType::Index3D getDim() const noexcept { return grid.getDim(); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(ssize_t x, ssize_t y, ssize_t z) const noexcept;
        [[nodiscard]] ScalarType getKinetic(ssize_t x, ssize_t y, ssize_t z) const noexcept;
    };

    template<class ScalarType>
    PWBaseWave<ScalarType>::PWBaseWave(ScalarType cutEnergy, LatticeMatrix repCell_) : repCell(std::move(repCell_)) {
        grid = KSpaceGrid<ComplexType>::makeGrid(cutEnergy, repCell);
    }

    template<class ScalarType>
    typename PWBaseWave<ScalarType>::ComplexType PWBaseWave<ScalarType>::operator()(Vector<ScalarType, 3> r) const {
        using Index3D = typename GridType::Index3D;
        ComplexType result = ComplexType::Zero();
        GridType::forKIndexInGrid(grid.getDim(), repCell, [this, &result, &r](Vector<ScalarType, 3> k, Index3D index) {
            const ScalarType phase = k * r;
            ScalarType s, c;
            sincos(phase, s, c);
            result += grid(index) * ComplexType(c, s);
        });
        return result;
    }

    template<class ScalarType>
    template<class VectorType>
    PWBaseWave<ScalarType>& PWBaseWave<ScalarType>::operator=(const RValueVector<VectorType>& newCoeffs) {
        grid.asVector() = newCoeffs;
        return *this;
    }

    template<class ScalarType>
    Vector<ScalarType, 3> PWBaseWave<ScalarType>::getWaveVector(ssize_t x, ssize_t y, ssize_t z) const noexcept {
        return repCell.row(0).asVector() * ScalarType(x) +
               repCell.row(1).asVector() * ScalarType(y) +
               repCell.row(2).asVector() * ScalarType(z);
    }

    template<class ScalarType>
    ScalarType PWBaseWave<ScalarType>::getKinetic(ssize_t x, ssize_t y, ssize_t z) const noexcept {
        constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
        return getWaveVector(x, y, z).squaredNorm() * ScalarType(factor);
    }
}
