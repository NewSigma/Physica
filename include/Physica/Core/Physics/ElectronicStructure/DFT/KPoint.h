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

#include "KPointBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"

namespace Physica::Core {
    template<class ScalarType, bool isSpinPolarized> class KPoint;

    template<class ScalarType>
    class KPoint<ScalarType, true> : public KPointBase<ScalarType> {
        using Base = KPointBase<ScalarType>;
        using typename Base::ComplexType;
        using typename Base::Vector3D;
        using typename Base::OccupacyVector;
        using EigInfo = EigenSolver<DenseMatrix<ComplexType>>;

        EigInfo eigUp;
        EigInfo eigDown;
        OccupacyVector occupacyUp;
        OccupacyVector occupacyDown;
    public:
        KPoint() = default;
        KPoint(Vector3D pos_, size_t plainWaveCount, ScalarType weight);
        KPoint(const KPoint&) = default;
        KPoint(KPoint&&) noexcept = default;
        ~KPoint() = default;
        /* Operators */
        KPoint& operator=(KPoint k) noexcept;
        /* Setters */
        void setEigInfo(EigInfo& newEigUp, EigInfo& newEigDown);
        /* Getters */
        [[nodiscard]] const EigInfo& getEigUp() const noexcept { return eigUp; }
        [[nodiscard]] const EigInfo& getEigDown() const noexcept { return eigDown; }
        [[nodiscard]] const OccupacyVector& getOccupacyUp() const noexcept { return occupacyUp; }
        [[nodiscard]] const OccupacyVector& getOccupacyDown() const noexcept { return occupacyDown; }
        /* Helpers */
        void swap(KPoint& kPoint) noexcept;
    };

    template<class ScalarType>
    KPoint<ScalarType, true>::KPoint(Vector3D pos_, size_t plainWaveCount, ScalarType weight)
            : Base(std::move(pos_), std::move(weight))
            , eigUp(plainWaveCount)
            , eigDown(plainWaveCount)
            , occupacyUp(plainWaveCount, 0)
            , occupacyDown(plainWaveCount, 0) {}

    template<class ScalarType>
    KPoint<ScalarType, true>& KPoint<ScalarType, true>::operator=(KPoint k) noexcept {
        swap(k);
        return *this;
    }

    template<class ScalarType>
    void KPoint<ScalarType, true>::setEigInfo(EigInfo& newEigUp, EigInfo& newEigDown) {
        std::swap(eigUp, newEigUp);
        std::swap(eigDown, newEigDown);
    }

    template<class ScalarType>
    void KPoint<ScalarType, true>::swap(KPoint& kPoint) noexcept {
        Base::swap(kPoint);
        std::swap(eigUp, kPoint.eigUp);
        std::swap(eigDown, kPoint.eigDown);
    }

    template<class ScalarType>
    class KPoint<ScalarType, false> : public KPointBase<ScalarType> {
        using Base = KPointBase<ScalarType>;
        using typename Base::ComplexType;
        using typename Base::Vector3D;
        using typename Base::OccupacyVector;
        using EigInfo = EigenSolver<DenseMatrix<ComplexType>>;

        EigInfo eig;
        OccupacyVector occupacy;
    public:
        KPoint() = default;
        KPoint(Vector3D pos_, size_t plainWaveCount, ScalarType weight);
        KPoint(const KPoint&) = default;
        KPoint(KPoint&&) noexcept = default;
        ~KPoint() = default;
        /* Operators */
        KPoint& operator=(KPoint k) noexcept;
        /* Setters */
        void setEigInfo(EigInfo& newEig);
        /* Getters */
        [[nodiscard]] const EigInfo& getEig() const noexcept { return eig; }
        [[nodiscard]] const OccupacyVector& getOccupacy() const noexcept { return occupacy; }
        /* Helpers */
        void swap(KPoint& kPoint) noexcept;
    };

    template<class ScalarType>
    KPoint<ScalarType, false>::KPoint(Vector3D pos_, size_t plainWaveCount, ScalarType weight)
            : Base(std::move(pos_), std::move(weight))
            , eig(plainWaveCount)
            , occupacy(plainWaveCount, 0) {}

    template<class ScalarType>
    KPoint<ScalarType, false>& KPoint<ScalarType, false>::operator=(KPoint k) noexcept {
        swap(k);
        return *this;
    }

    template<class ScalarType>
    void KPoint<ScalarType, false>::setEigInfo(EigInfo& newEig) {
        std::swap(eig, newEig);
    }

    template<class ScalarType>
    void KPoint<ScalarType, false>::swap(KPoint& kPoint) noexcept {
        Base::swap(kPoint);
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
