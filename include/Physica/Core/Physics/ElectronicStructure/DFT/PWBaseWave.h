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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Grid3D.h"

namespace Physica::Core {
    /**
     * \class PWBaseWave is plain wave based wave function
     */
    template<class ScalarType>
    class PWBaseWave {
        using ComplexType = ComplexScalar<ScalarType>;
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using SignedGrid = Grid3D<ComplexType, true>;
        using Dim = typename SignedGrid::Dim;

        SignedGrid grid;
    public:
        PWBaseWave(ScalarType cutEnergy, LatticeMatrix reciprocalLattice);
        /* Operators */
        [[nodiscard]] ComplexType operator()(Vector<ScalarType, 3> r) const;
        template<class VectorType>
        PWBaseWave& operator=(const RValueVector<VectorType>& newCoeffs);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return grid.getSize(); }
        [[nodiscard]] size_t getDimX() const noexcept { return grid.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return grid.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return grid.getDimZ(); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(Dim dim) const noexcept;
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(size_t index) const noexcept { return getWaveVector(indexToDim(index)); }
        [[nodiscard]] ScalarType getKinetic(Dim dim) const noexcept;
        [[nodiscard]] size_t dimToIndex(ssize_t x, ssize_t y, ssize_t z) const noexcept { return grid.dimToIndex(x, y, z); }
        [[nodiscard]] Dim indexToDim(size_t index) const noexcept { return grid.indexToDim(index); }
    };

    template<class ScalarType>
    PWBaseWave<ScalarType>::PWBaseWave(ScalarType cutEnergy, LatticeMatrix reciprocalLattice)
            : grid(SignedGrid::gridFromCutEnergy(cutEnergy, reciprocalLattice)) {}

    template<class ScalarType>
    typename PWBaseWave<ScalarType>::ComplexType PWBaseWave<ScalarType>::operator()(Vector<ScalarType, 3> r) const {
        const size_t length = grid.getSize();
        ComplexType result = ComplexType::Zero();
        for (size_t i = 0; i < length; ++i) {
            const ScalarType phase = getWaveVector(indexToDim(i)) * r;
            const auto c = cos(phase);
            const auto s = sqrt(ScalarType::One() - square(c)); //Sign is not necessary to wave functions
            result += grid[i] * ComplexType(c, s); //Using platform specific sincos() may gain performance improvement
        }
        return result;
    }

    template<class ScalarType>
    template<class VectorType>
    PWBaseWave<ScalarType>& PWBaseWave<ScalarType>::operator=(const RValueVector<VectorType>& newCoeffs) {
        grid.asVector() = newCoeffs;
        return *this;
    }

    template<class ScalarType>
    Vector<ScalarType, 3> PWBaseWave<ScalarType>::getWaveVector(Dim dim) const noexcept {
        auto[n1, n2, n3] = dim;
        const auto& lattice = grid.getLattice();
        return lattice.row(0).asVector() * ScalarType(n1) +
               lattice.row(1).asVector() * ScalarType(n2) +
               lattice.row(2).asVector() * ScalarType(n3);
    }

    template<class ScalarType>
    ScalarType PWBaseWave<ScalarType>::getKinetic(Dim dim) const noexcept {
        constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
        return getWaveVector(dim).squaredNorm() * ScalarType(factor);
    }
}
