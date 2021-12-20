/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType>
    class WaveFunction {
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using SignedGrid = Grid3D<ScalarType, true>;
        using Dim = typename SignedGrid::Dim;

        SignedGrid grid;
    public:
        WaveFunction(ScalarType cutEnergy, LatticeMatrix mat);
        /* Operators */
        [[nodiscard]] ComplexScalar<ScalarType> operator()(Vector<ScalarType, 3> r) const;
        template<class VectorType>
        WaveFunction& operator=(const RValueVector<VectorType>& newCoeffs);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return grid.getSize(); }
        [[nodiscard]] size_t getDimX() const noexcept { return grid.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return grid.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return grid.getDimZ(); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(Dim dim) const noexcept;
        [[nodiscard]] ScalarType getKinetic(Dim dim) const noexcept;
        [[nodiscard]] size_t dimToIndex(ssize_t x, ssize_t y, ssize_t z) const noexcept { return grid.dimToIndex(x, y, z); }
        [[nodiscard]] Dim indexToDim(size_t index) const noexcept { return grid.indexToDim(index); }
    };

    template<class ScalarType>
    WaveFunction<ScalarType>::WaveFunction(ScalarType cutEnergy, LatticeMatrix mat) : grid() {
        const ScalarType factor = ScalarType(2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck);
        const ScalarType maxMoment = sqrt(factor * cutEnergy.getTrivial());
        const auto dimX = size_t((maxMoment / mat.row(0).norm()).getTrivial());
        const auto dimY = size_t((maxMoment / mat.row(1).norm()).getTrivial());
        const auto dimZ = size_t((maxMoment / mat.row(2).norm()).getTrivial());
        
        grid = SignedGrid(mat, dimX, dimY, dimZ);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> WaveFunction<ScalarType>::operator()(Vector<ScalarType, 3> r) const {
        const size_t length = grid.getSize();
        ComplexScalar<ScalarType> result = ComplexScalar<ScalarType>::Zero();
        for (size_t i = 0; i < length; ++i) {
            const ScalarType phase = getWaveVector(indexToDim(i)) * r;
            result += grid[i] * ComplexScalar<ScalarType>(cos(phase), sin(phase));
        }
        return result;
    }

    template<class ScalarType>
    template<class VectorType>
    WaveFunction<ScalarType>& WaveFunction<ScalarType>::operator=(const RValueVector<VectorType>& newCoeffs) {
        grid.asVector() = newCoeffs;
        return *this;
    }

    template<class ScalarType>
    Vector<ScalarType, 3> WaveFunction<ScalarType>::getWaveVector(Dim dim) const noexcept {
        auto[n1, n2, n3] = dim;
        const auto& lattice = grid.getLattice();
        return lattice.row(0).asVector() * ScalarType(n1) +
               lattice.row(1).asVector() * ScalarType(n2) +
               lattice.row(2).asVector() * ScalarType(n3);
    }

    template<class ScalarType>
    ScalarType WaveFunction<ScalarType>::getKinetic(Dim dim) const noexcept {
        constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
        return getWaveVector(dim).squaredNorm() * ScalarType(factor);
    }
}
