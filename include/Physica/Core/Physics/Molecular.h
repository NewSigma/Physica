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

#include <fstream>
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Math/Geometry/Point.h"
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core::Physics {
    template<class ScalarType>
    class Molecular {
        using PointType = Point<3, ScalarType>;
        Array<PointType> atoms;
        Array<unsigned char> atomicNumbers;
    public:
        Molecular();
        ~Molecular() = default;
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const Molecular& m);
        /* Getters */
        [[nodiscard]] size_t atomCount() const noexcept { return atoms.getLength(); }
        [[nodiscard]] PointType getAtom(size_t i) const { return atoms[i]; }
        [[nodiscard]] unsigned char getAtomicNumber(size_t i) const { return atomicNumber[i]; }
        [[nodiscard]] ScalarType getPairDist(size_t i, size_t j) const;
        [[nodiscard]] ScalarType getTripleAngle(size_t i, size_t j, size_t k) const;
        [[nodiscard]] ScalarType getOutOfPlaneAngle(size_t i, size_t j, size_t k, size_t l) const;
        [[nodiscard]] ScalarType getDihedralAngle(size_t i, size_t j, size_t k, size_t l) const;
        [[nodiscard]] PointType getMassCenter() const;
    };

    template<class ScalarType>
    std::ostream& operator<<(std::ostream& os, const Molecular& m) {
        os << atoms.getLength() << '\n';
        for (auto& atom : atoms)
            os << atom << '\n';
        return os;
    }

    template<class ScalarType>
    ScalarType Molecular<ScalarType>::getPairDist(size_t i, size_t j) const {
        return atoms[i].dist(atoms[j]);
    }
    /**
     * \returns The angle between vector ji and vector jk
     */
    template<class ScalarType>
    ScalarType Molecular<ScalarType>::getTripleAngle(size_t i, size_t j, size_t k) {
        using VectorType = typename PointType::VectorType;
        const VectorType vector_ji = atoms[i].getVector() - atoms[j].getVector();
        const VectorType vector_jk = atoms[k].getVector() - atoms[j].getVector();
        return arccos(vector_ji * vector_jk / (vector_ji.norm() * vector_jk.norm()));
    }
    /**
     * \returns Angle between plane ijk and line kl
     */
    template<class ScalarType>
    ScalarType Molecular<ScalarType>::getOutOfPlaneAngle(size_t i, size_t j, size_t k, size_t l) const {
        using VectorType = typename PointType::VectorType;
        const VectorType vector_ki = atoms[i].getVector() - atoms[k].getVector();
        const VectorType vector_kj = atoms[j].getVector() - atoms[k].getVector();
        const VectorType cross = crossProduct(vector_ki, vector_kj);
        const VectorType vector_kl = atoms[l].getVector() - atoms[k].getVector();
        return cross * vector_kl / (cross.norm() * vector_kl.norm());
    }
    /**
     * \returns Dihedral angle between plain ijk and plain jkl
     */
    template<class ScalarType>
    ScalarType Molecular<ScalarType>::getDihedralAngle(size_t i, size_t j, size_t k, size_t l) const {
        using VectorType = typename PointType::VectorType;
        const VectorType vector_ki = atoms[i].getVector() - atoms[k].getVector();
        const VectorType vector_kj = atoms[j].getVector() - atoms[k].getVector();
        const VectorType cross1 = crossProduct(vector_ki, vector_kj);
        const VectorType vector_kl = atoms[l].getVector() - atoms[k].getVector();
        const VectorType cross2 = crossProduct(vector_kl, vector_kj);
        return cross1 * cross2 / (cross1.norm() * cross2.norm());
    }

    template<class ScalarType>
    Molecular<ScalarType>::PointType Molecular<ScalarType>::getMassCenter() const {
        double totalMass = 0;
        const size_t length = atoms.getLength();
        VectorType result = VectorType::Zeros(length);
        for (size_t i = 0; i < length; ++i) {
            double atomMass = PhyConst::relativeAtomMass[atomicNumbers[i]];
            totalMass += atomMass;
            result += ScalarType(atomMass) * atoms[i].getVector();
        }
        result *= reciprocal(ScalarType(totalMass));
        return PointType(std::move(result));
    }
}
