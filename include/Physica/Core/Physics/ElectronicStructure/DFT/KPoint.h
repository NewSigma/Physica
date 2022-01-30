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
    template<class ScalarType, bool isSpinPolarized> class KPoint;

    template<class ScalarType>
    class KPoint<ScalarType, true> {
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using EigInfo = EigenSolver<DenseMatrix<ComplexType>>;

        Vector3D pos;
        EigInfo eigUp;
        EigInfo eigDown;
    public:
        KPoint() = default;
        KPoint(Vector3D pos_, size_t plainWaveCount);
        KPoint(const KPoint&) = default;
        KPoint(KPoint&&) = default;
        ~KPoint() = default;
        /* Operators */
        KPoint& operator=(const KPoint&) = default;
        KPoint& operator=(KPoint&&) = default;
        /* Setters */
        void setEigInfo(EigInfo& newEigUp, EigInfo& newEigDown);
        /* Getters */
        [[nodiscard]] const Vector3D& getPos() { return pos; }
        [[nodiscard]] const EigInfo& getEigUp() const noexcept { return eigUp; }
        [[nodiscard]] const EigInfo& getEigDown() const noexcept { return eigDown; }
        /* Helpers */
        void swap(KPoint& kPoint) noexcept;
    };

    template<class ScalarType>
    KPoint<ScalarType, true>::KPoint(Vector3D pos_, size_t plainWaveCount)
            : pos(std::move(pos_)), eigUp(plainWaveCount), eigDown(plainWaveCount) {}

    template<class ScalarType>
    void KPoint<ScalarType, true>::setEigInfo(EigInfo& newEigUp, EigInfo& newEigDown) {
        std::swap(eigUp, newEigUp);
        std::swap(eigDown, newEigDown);
    }

    template<class ScalarType>
    void KPoint<ScalarType, true>::swap(KPoint& kPoint) noexcept {
        std::swap(pos, kPoint.pos);
        std::swap(eigUp, kPoint.eigUp);
        std::swap(eigDown, kPoint.eigDown);
    }

    template<class ScalarType>
    class KPoint<ScalarType, false> {
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using EigInfo = EigenSolver<DenseMatrix<ComplexType>>;

        Vector3D pos;
        EigInfo eig;
    public:
        KPoint() = default;
        KPoint(Vector3D pos_, size_t plainWaveCount);
        KPoint(const KPoint&) = default;
        KPoint(KPoint&&) = default;
        ~KPoint() = default;
        /* Operators */
        KPoint& operator=(const KPoint&) = default;
        KPoint& operator=(KPoint&&) = default;
        /* Setters */
        void setEigInfo(EigInfo& newEig);
        /* Getters */
        [[nodiscard]] const Vector3D& getPos() { return pos; }
        [[nodiscard]] const EigInfo& getEig() const noexcept { return eig; }
        /* Helpers */
        void swap(KPoint& kPoint) noexcept;
    };

    template<class ScalarType>
    KPoint<ScalarType, false>::KPoint(Vector3D pos_, size_t plainWaveCount)
            : pos(std::move(pos_)), eig(plainWaveCount) {}

    template<class ScalarType>
    void KPoint<ScalarType, false>::setEigInfo(EigInfo& newEig) {
        std::swap(eig, newEig);
    }

    template<class ScalarType>
    void KPoint<ScalarType, false>::swap(KPoint& kPoint) noexcept {
        std::swap(pos, kPoint.pos);
        std::swap(eig, kPoint.eig);
    }
}

namespace std {
    template<class ScalarType, bool isSpinPolarized>
    inline void swap(Physica::Core::KPoint<ScalarType, isSpinPolarized>& kPoint1,
                     Physica::Core::KPoint<ScalarType, isSpinPolarized>& kPoint2) noexcept {
        kPoint1.swap(kPoint2);
    }
}
