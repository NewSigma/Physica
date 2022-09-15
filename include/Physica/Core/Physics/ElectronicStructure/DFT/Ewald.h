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

#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Grid3D.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] VASP (www.vasp.at)
     */
    template<class OutScalarType, class PosScalarType = OutScalarType>
    class Ewald {
        constexpr static unsigned int Dim = 3;
        using LatticeMatrix = typename PeriodicCell<PosScalarType, Dim>::LatticeMatrix;
        using PositionMatrix = typename PeriodicCell<PosScalarType, Dim>::PositionMatrix;
        using Vector3D = Vector<PosScalarType, Dim>;

        LatticeMatrix lattice;
        ReciprocalCell<PosScalarType> repCell;
        Vector<OutScalarType> charges;
        PosScalarType inv_volume;
        PosScalarType averageCellSize;
        OutScalarType integralLimit;
        Utils::Array<int, Dim> realSumRange;
        Utils::Array<int, Dim> repSumRange;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice_, Vector<OutScalarType> charges_);
        Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<OutScalarType> charges_);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald ewald) noexcept;
        [[nodiscard]] OutScalarType operator()(const PositionMatrix& pos) const;
        /* Operations */
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
    private:
        OutScalarType realSum(const Vector3D& deltaPos) const;
        OutScalarType reciprocalSum(const Vector3D& deltaPos) const;
        static Utils::Array<int, Dim> getSumRange(const LatticeMatrix& lattice, PosScalarType factor);
    };

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, Vector<OutScalarType> charges_)
            : Ewald(lattice_, ReciprocalCell(lattice_), std::move(charges_)) {}

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<OutScalarType> charges_)
            : lattice(std::move(lattice_))
            , repCell(std::move(repCell_))
            , charges(std::move(charges_)) {
        const PosScalarType volume = PeriodicCell<PosScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(volume);
        //The following param chosen is referenced from VASP
        averageCellSize = cbrt(volume);
        integralLimit = sqrt(OutScalarType(M_PI)) / averageCellSize;
        realSumRange = getSumRange(repCell.getLattice(), OutScalarType(2 / M_PI) / integralLimit);
        repSumRange = getSumRange(lattice, OutScalarType(4 / M_PI) * integralLimit);
    }

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>& Ewald<OutScalarType, PosScalarType>::operator=(Ewald ewald) noexcept {
        swap(ewald);
        return *this;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::operator()(const PositionMatrix& pos) const {
        using ScalarType = OutScalarType;
        ScalarType result = ScalarType::Zero();
        ScalarType totalCharge = 0;
        ScalarType totalSquaredCharge = 0;
        for (size_t i = 0; i < getNumParticle(); ++i) {
            for (size_t j = 0; j < getNumParticle(); ++j) { //Optimize: possible to loop from ion2 = ion1
                const Vector<PosScalarType, Dim> deltaPos = pos.row(i).asVector() - pos.row(j);
                ScalarType sum = realSum(deltaPos); //Optimize: VASP puts this loop outside, consider its performance
                sum += ScalarType(4 * M_PI) * reciprocalSum(deltaPos) * inv_volume;

                const ScalarType dotCharge = charges[i] * charges[j];
                result += sum * dotCharge;
            }
            const ScalarType charge = charges[i];
            totalCharge += charge;
            totalSquaredCharge += square(charge);
        }
        result *= ScalarType(0.5);
        result -= totalSquaredCharge * integralLimit / sqrt(ScalarType(M_PI));
        result -= square(totalCharge) * ScalarType(M_PI) / (ScalarType::Two() * square(integralLimit)) * inv_volume;
        return result;
    }

    template<class OutScalarType, class PosScalarType>
    void Ewald<OutScalarType, PosScalarType>::swap(Ewald& ewald) noexcept {
        lattice.swap(ewald.lattice);
        repCell.swap(ewald.repCell);
        charges.swap(ewald.charges);
        inv_volume.swap(ewald.inv_volume);
        averageCellSize.swap(ewald.averageCellSize);
        integralLimit.swap(ewald.integralLimit);
        realSumRange.swap(ewald.realSumRange);
        repSumRange.swap(ewald.repSumRange);
    }

    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::realSum(const Vector3D& deltaPos) const {
        using ScalarType = OutScalarType;
        ScalarType sum = ScalarType::Zero();

        const bool isTooClose = deltaPos.norm() <= averageCellSize * std::numeric_limits<ScalarType>::epsilon();
        if (isTooClose) {
            for (int i = -realSumRange[0]; i <= realSumRange[0]; ++i) {
                for (int j = -realSumRange[1]; j <= realSumRange[1]; ++j) {
                    for (int k = -realSumRange[2]; k <= realSumRange[2]; ++k) {
                        if (i == 0 && j == 0 && k == 0)
                            continue;
                        const Vector3D transVector = Vector3D(lattice.row(0)) * ScalarType(i)
                                                   + Vector3D(lattice.row(1)) * ScalarType(j)
                                                   + Vector3D(lattice.row(2)) * ScalarType(k);
                        const ScalarType norm = transVector.norm();
                        sum += erfc(integralLimit * norm) / norm; //Optimize: VASP uses searching table method
                    }
                }
            }
        }
        else {
            for (int i = -realSumRange[0]; i <= realSumRange[0]; ++i) {
                for (int j = -realSumRange[1]; j <= realSumRange[1]; ++j) {
                    for (int k = -realSumRange[2]; k <= realSumRange[2]; ++k) {
                        const Vector3D transVector = Vector3D(lattice.row(0)) * ScalarType(i)
                                                   + Vector3D(lattice.row(1)) * ScalarType(j)
                                                   + Vector3D(lattice.row(2)) * ScalarType(k);
                        const ScalarType norm = (deltaPos - transVector).norm();
                        sum += erfc(integralLimit * norm) / norm;
                    }
                }
            }
        }
        return sum;
    }

    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::reciprocalSum(const Vector3D& deltaPos) const {
        using ScalarType = OutScalarType;
        const auto& repLattice = repCell.getLattice();
        ScalarType sum = ScalarType::Zero();
        for (int i = -repSumRange[0]; i <= repSumRange[0]; ++i) {
            for (int j = -repSumRange[1]; j <= repSumRange[1]; ++j) {
                for (int k = -repSumRange[2]; k <= repSumRange[2]; ++k) {
                    if (i == 0 && j == 0 && k == 0)
                        continue;
                    const Vector3D repVector = Vector3D(repLattice.row(0)) * ScalarType(i)
                                             + Vector3D(repLattice.row(1)) * ScalarType(j)
                                             + Vector3D(repLattice.row(2)) * ScalarType(k);
                    const ScalarType squaredNorm = repVector.squaredNorm();
                    const ScalarType dot = repVector * deltaPos;
                    sum += cos(dot) / (squaredNorm * exp(squaredNorm / square(ScalarType::Two() * integralLimit)));
                }
            }
        }
        return sum;
    }

    template<class OutScalarType, class PosScalarType>
    Utils::Array<int, Ewald<OutScalarType, PosScalarType>::Dim>
    Ewald<OutScalarType, PosScalarType>::getSumRange(const LatticeMatrix& lattice, PosScalarType factor) {
        constexpr double roundFactor = 1 - std::numeric_limits<double>::epsilon();
        static_assert(roundFactor < 1);

        const int dim1 = int((factor * lattice.row(0).norm() + roundFactor).getTrivial());
        const int dim2 = int((factor * lattice.row(1).norm() + roundFactor).getTrivial());
        const int dim3 = int((factor * lattice.row(2).norm() + roundFactor).getTrivial());
        return {dim1, dim2, dim3};
    }
}
