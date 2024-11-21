/*
 * Copyright 2021-2024 Weibo He.
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

#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Math/Geometry/Point.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica::Core::Physics {
    template<Scalar T>
    class Molecular {
        using PointType = Point<3, T>;
        Array<PointType> atoms;
        Array<uint8_t> atomicNumbers;
    public:
        Molecular(size_t atomCount);
        ~Molecular() = default;
        /* Getters */
        [[nodiscard]] Array<PointType>& getAtoms() noexcept { return atoms; }
        [[nodiscard]] Array<uint8_t>& getAtomicNumbers() noexcept { return atomicNumbers; }
        [[nodiscard]] size_t getAtomCount() const noexcept { return atoms.getLength(); }
        [[nodiscard]] PointType getAtom(size_t i) const { return atoms[i]; }
        [[nodiscard]] uint8_t getAtomicNumber(size_t i) const { return atomicNumbers[i]; }
        [[nodiscard]] T getPairDist(size_t i, size_t j) const;
        [[nodiscard]] T getTripleAngle(size_t i, size_t j, size_t k) const;
        [[nodiscard]] T getOutOfPlaneAngle(size_t i, size_t j, size_t k, size_t l) const;
        [[nodiscard]] T getDihedralAngle(size_t i, size_t j, size_t k, size_t l) const;
        [[nodiscard]] PointType getMassCenter() const;
        [[nodiscard]] T getNuclearRepulsionEnergy() const;
    };

    template<Scalar T>
    Molecular<T>::Molecular(size_t atomCount) : atoms(atomCount), atomicNumbers(atomCount) {}

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const Molecular<T>& m) {
        os << m.atomCount() << '\n';
        const auto& atoms = m.getAtoms();
        for (const auto& atom : atoms)
            os << atom << '\n';
        return os;
    }

    template<Scalar T>
    T Molecular<T>::getPairDist(size_t i, size_t j) const {
        return atoms[i].dist(atoms[j]);
    }
    /**
     * \returns The angle between vector ji and vector jk
     */
    template<Scalar T>
    T Molecular<T>::getTripleAngle(size_t i, size_t j, size_t k) const {
        using VectorType = PointType::VectorType;
        const VectorType vector_ji = atoms[i].getVector() - atoms[j].getVector();
        const VectorType vector_jk = atoms[k].getVector() - atoms[j].getVector();
        return arccos(vector_ji * vector_jk / (vector_ji.norm() * vector_jk.norm()));
    }
    /**
     * \returns Angle between plane ijk and line kl
     */
    template<Scalar T>
    T Molecular<T>::getOutOfPlaneAngle(size_t i, size_t j, size_t k, size_t l) const {
        using VectorType = PointType::VectorType;
        const VectorType vector_ki = atoms[i].getVector() - atoms[k].getVector();
        const VectorType vector_kj = atoms[j].getVector() - atoms[k].getVector();
        const VectorType cross = crossProduct(vector_ki, vector_kj);
        const VectorType vector_kl = atoms[l].getVector() - atoms[k].getVector();
        return cross * vector_kl / (cross.norm() * vector_kl.norm());
    }
    /**
     * \returns Dihedral angle between plain ijk and plain jkl
     */
    template<Scalar T>
    T Molecular<T>::getDihedralAngle(size_t i, size_t j, size_t k, size_t l) const {
        using VectorType = PointType::VectorType;
        const VectorType vector_ki = atoms[i].getVector() - atoms[k].getVector();
        const VectorType vector_kj = atoms[j].getVector() - atoms[k].getVector();
        const VectorType cross1 = crossProduct(vector_ki, vector_kj);
        const VectorType vector_kl = atoms[l].getVector() - atoms[k].getVector();
        const VectorType cross2 = crossProduct(vector_kl, vector_kj);
        return cross1 * cross2 / (cross1.norm() * cross2.norm());
    }

    template<Scalar T>
    Molecular<T>::PointType Molecular<T>::getMassCenter() const {
        using VectorType = PointType::VectorType;
        double totalMass = 0;
        const size_t length = atoms.getLength();
        auto result = VectorType::zeros(length);
        for (size_t i = 0; i < length; ++i) {
            double atomMass = PhyConst<SI>::relativeAtomMass[atomicNumbers[i]];
            totalMass += atomMass;
            result += T(atomMass) * atoms[i].getVector();
        }
        result *= reciprocal(T(totalMass));
        return PointType(std::move(result));
    }
    /**
     * Using atom unit
     */
    template<Scalar T>
    T Molecular<T>::getNuclearRepulsionEnergy() const {
        T result = T(0);
        const size_t atomCount = getAtomCount();
        for (size_t i = 0; i < atomCount - 1; ++i) {
            const T atomicNumber1 = T(getAtomicNumber(i));
            for (size_t j = i + 1; j < atomCount; ++j) {
                const T atomicNumber2 = T(getAtomicNumber(j));
                const T dist = getPairDist(i, j);
                result += atomicNumber1 * atomicNumber2 / dist;
            }
        }
        return result;
    }
}
