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

namespace Physica::Core {
    template<class ScalarType>
    class WaveFunction {
        using LatticeMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, 3, 3>;
        using Dim = std::tuple<ssize_t, ssize_t, ssize_t>;

        LatticeMatrix reciprocalLattice;
        /**
         * Coeff makes a grid of size [-baseDimX, baseDimX] * [-baseDimY, baseDimY] * [-baseDimZ, baseDimZ]
         */
        Vector<ScalarType> baseCoeffs;
        size_t baseDimX;
        size_t baseDimY;
        size_t baseDimZ;
    public:
        WaveFunction(ScalarType cutEnergy, LatticeMatrix mat);
        /* Operators */
        [[nodiscard]] ComplexScalar<ScalarType> operator()(Vector<ScalarType, 3> r) const;
        template<class VectorType>
        WaveFunction& operator=(const RValueVector<VectorType>& newCoeffs);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return baseCoeffs.getLength(); }
        [[nodiscard]] size_t getDimX() const noexcept { return baseDimX; }
        [[nodiscard]] size_t getDimY() const noexcept { return baseDimY; }
        [[nodiscard]] size_t getDimZ() const noexcept { return baseDimZ; }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(Dim dim) const noexcept;
        [[nodiscard]] ScalarType getKinetic(Dim dim) const noexcept;
        [[nodiscard]] size_t dimToIndex(ssize_t x, ssize_t y, ssize_t z) const noexcept;
        [[nodiscard]] Dim indexToDim(size_t index) const noexcept;
    };

    template<class ScalarType>
    WaveFunction<ScalarType>::WaveFunction(ScalarType cutEnergy, LatticeMatrix mat) : reciprocalLattice(std::move(mat)) {
        using namespace Physica::Core::Physics;
        const ScalarType factor = ScalarType(2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck);
        const ScalarType maxMoment = sqrt(factor * cutEnergy.getTrivial());
        baseDimX = size_t((maxMoment / reciprocalLattice.row(0).norm()).getTrivial());
        baseDimY = size_t((maxMoment / reciprocalLattice.row(1).norm()).getTrivial());
        baseDimZ = size_t((maxMoment / reciprocalLattice.row(2).norm()).getTrivial());
        
        const size_t baseSize = (baseDimX * 2 + 1) * (baseDimY * 2 + 1) * (baseDimZ * 2 + 1);
        baseCoeffs.resize(baseSize);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> WaveFunction<ScalarType>::operator()(Vector<ScalarType, 3> r) const {
        const size_t length = baseCoeffs.getLength();
        ComplexScalar<ScalarType> result = ComplexScalar<ScalarType>::Zero();
        for (size_t i = 0; i < length; ++i) {
            const ScalarType phase = getWaveVector(indexToDim(i)) * r;
            result += baseCoeffs[i] * ComplexScalar<ScalarType>(cos(phase), sin(phase));
        }
        return result;
    }

    template<class ScalarType>
    template<class VectorType>
    WaveFunction<ScalarType>& WaveFunction<ScalarType>::operator=(const RValueVector<VectorType>& newCoeffs) {
        baseCoeffs = newCoeffs;
        return *this;
    }

    template<class ScalarType>
    Vector<ScalarType, 3> WaveFunction<ScalarType>::getWaveVector(Dim dim) const noexcept {
        auto[n1, n2, n3] = dim;
        return reciprocalLattice.row(0).asVector() * ScalarType(n1) +
               reciprocalLattice.row(1).asVector() * ScalarType(n2) +
               reciprocalLattice.row(2).asVector() * ScalarType(n3);
    }

    template<class ScalarType>
    ScalarType WaveFunction<ScalarType>::getKinetic(Dim dim) const noexcept {
        using namespace Physica::Core::Physics;
        constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
        return getWaveVector(dim).squaredNorm() * ScalarType(factor);
    }

    template<class ScalarType>
    size_t WaveFunction<ScalarType>::dimToIndex(ssize_t x, ssize_t y, ssize_t z) const noexcept {
        const size_t sizeX = 2 * baseDimX + 1;
        const size_t sizeXY = sizeX * (2 * baseDimY + 1);
        
        const size_t normalized_x = x + static_cast<ssize_t>(baseDimX);
        const size_t normalized_y = y + static_cast<ssize_t>(baseDimY);
        const size_t normalized_z = z + static_cast<ssize_t>(baseDimZ);
        return normalized_z * sizeXY + normalized_y * sizeX + normalized_x;
    }

    template<class ScalarType>
    typename WaveFunction<ScalarType>::Dim WaveFunction<ScalarType>::indexToDim(size_t index) const noexcept {
        const size_t sizeX = 2 * baseDimX + 1;
        const size_t sizeXY = sizeX * (2 * baseDimY + 1);

        const size_t normalized_z = index / sizeXY;
        index -= sizeXY * normalized_z;
        const size_t normalized_y = index / sizeX;
        index -= sizeX * normalized_y;
        const size_t normalized_x = index;

        return {static_cast<ssize_t>(normalized_x) - static_cast<ssize_t>(baseDimX),
                static_cast<ssize_t>(normalized_y) - static_cast<ssize_t>(baseDimY),
                static_cast<ssize_t>(normalized_z) - static_cast<ssize_t>(baseDimZ)};
    }
}
