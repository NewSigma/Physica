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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"

namespace Physica::Core {
    template<class ScalarType, bool isSpinPolarized> class KPointInfo;

    template<class ScalarType>
    class KPointInfo<ScalarType, true> {
        using ComplexType = ComplexScalar<ScalarType>;
        using EigInfo = EigenSolver<DenseMatrix<ComplexType>>;

        EigInfo eigUp;
        EigInfo eigDown;
    public:
        KPointInfo(size_t plainWaveCount);
        /* Setters */
        void setEigInfo(EigInfo& newEigUp, EigInfo& newEigDown);
        /* Getters */
        const EigInfo& getEigUp() const noexcept { return eigUp; }
        const EigInfo& getEigDown() const noexcept { return eigDown; }
    };

    template<class ScalarType>
    KPointInfo<ScalarType, true>::KPointInfo(size_t plainWaveCount) : eigUp(plainWaveCount), eigDown(plainWaveCount) {}

    template<class ScalarType>
    void KPointInfo<ScalarType, true>::setEigInfo(EigInfo& newEigUp, EigInfo& newEigDown) {
        std::swap(eigUp, newEigUp);
        std::swap(eigDown, newEigDown);
    }

    template<class ScalarType>
    class KPointInfo<ScalarType, false> {
        using ComplexType = ComplexScalar<ScalarType>;
        using EigInfo = EigenSolver<DenseMatrix<ComplexType>>;

        EigInfo eig;
    public:
        KPointInfo(size_t plainWaveCount);
        /* Setters */
        void setEigInfo(EigInfo& newEig);
        /* Getters */
        const EigInfo& getEig() const noexcept { return eig; }
    };

    template<class ScalarType>
    KPointInfo<ScalarType, false>::KPointInfo(size_t plainWaveCount) : eig(plainWaveCount) {}

    template<class ScalarType>
    void KPointInfo<ScalarType, false>::setEigInfo(EigInfo& newEig) {
        std::swap(eig, newEig);
    }
}
