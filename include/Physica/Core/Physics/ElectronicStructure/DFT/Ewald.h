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
    template<class ScalarType>
    class Ewald {
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using RepCellType = ReciprocalCell<typename LatticeMatrix::ScalarType>;
        using UnsignedGrid = Grid3D<ScalarType, false>;
        using Vector3D = Vector<ScalarType, 3>;
    public:
        static ScalarType energyIonIon(CrystalCell cell, const RepCellType& repCell, const Utils::Array<int16_t>& charges);
        static ScalarType energyIonIon(const MDCell& cell, const RepCellType& repCell, const Vector<ScalarType>& charges);
        template<class ScalarType_>
        static ScalarType energyIonIon(const PeriodicCell<ScalarType_, 3>& cell,
                                       const RepCellType& repCell,
                                       const Vector<ScalarType>& charges);
    private:
        [[nodiscard]] static std::tuple<int, int, int> getSumDimention(const LatticeMatrix& latt, ScalarType factor);
        [[nodiscard]] static ScalarType realSum(const LatticeMatrix& cell,
                                                ScalarType integralLimit,
                                                std::tuple<int, int, int> dim,
                                                const Vector<ScalarType, 3>& deltaPos,
                                                const ScalarType& averageCellSize);
        [[nodiscard]] static ScalarType reciprocalSum(const LatticeMatrix& cell,
                                                      ScalarType integralLimit,
                                                      std::tuple<int, int, int> dim,
                                                      const Vector<ScalarType, 3>& deltaPos);
        [[nodiscard]] static bool isTooClose(ScalarType deltaR, ScalarType averageCellSize);
    };

    template<class ScalarType>
    ScalarType Ewald<ScalarType>::energyIonIon(CrystalCell cell, const RepCellType& repCell, const Utils::Array<int16_t>& charges) {
        if (cell.getType() == CrystalCell::Type::Direct)
            cell.toCartesian();
        Vector<ScalarType> temp(charges.getLength());
        for (size_t i = 0; i < charges.getLength(); ++i)
            temp[i] = ScalarType(charges[i]);
        return energyIonIon(cell, repCell, temp);
    }

    template<class ScalarType>
    ScalarType Ewald<ScalarType>::energyIonIon(const MDCell& cell, const RepCellType& repCell, const Vector<ScalarType>& charges) {
        return energyIonIon(cell, repCell, charges);
    }
    /**
     * \param cell must be in cartesian convension
     */
    template<class ScalarType>
    template<class ScalarType_>
    ScalarType Ewald<ScalarType>::energyIonIon(const PeriodicCell<ScalarType_, 3>& cell,
                                               const RepCellType& repCell,
                                               const Vector<ScalarType>& charges) {
        const size_t ionCount = cell.getNumParticle();
        const ScalarType inv_volume = reciprocal(cell.getVolume());
        //The following param chosen is referenced from VASP
        const ScalarType averageCellSize = cbrt(ScalarType(cell.getVolume()));
        const ScalarType integralLimit = sqrt(ScalarType(M_PI)) / averageCellSize;
        const auto realSumDim = getSumDimention(repCell.getLattice(), ScalarType(2 / M_PI) / integralLimit);
        const auto repSumDim = getSumDimention(cell.getLattice(), ScalarType(4 / M_PI) * integralLimit);

        ScalarType result = ScalarType::Zero();
        ScalarType totalCharge = 0;
        ScalarType totalSquaredCharge = 0;
        for (size_t ion1 = 0; ion1 < ionCount; ++ion1) {
            for (size_t ion2 = 0; ion2 < ionCount; ++ion2) { //Optimize: possible to loop from ion2 = ion1
                const Vector3D deltaPos = cell.getPos().row(ion1).asVector() - cell.getPos().row(ion2);
                ScalarType sum = realSum(cell.getLattice(), integralLimit, realSumDim, deltaPos, averageCellSize); //Optimize: VASP puts this loop outside, consider its performance
                sum += ScalarType(4 * M_PI) * reciprocalSum(repCell.getLattice(), integralLimit, repSumDim, deltaPos) * inv_volume;

                const ScalarType dotCharge = charges[ion1] * charges[ion2];
                result += sum * dotCharge;
            }
            const ScalarType charge = charges[ion1];
            totalCharge += charge;
            totalSquaredCharge += square(charge);
        }
        result *= ScalarType(0.5);
        result -= totalSquaredCharge * integralLimit / sqrt(ScalarType(M_PI));
        result -= square(totalCharge) * ScalarType(M_PI) / (ScalarType::Two() * square(integralLimit)) * inv_volume;
        return result;
    }

    template<class ScalarType>
    std::tuple<int, int, int> Ewald<ScalarType>::getSumDimention(const LatticeMatrix& latt, ScalarType factor) {
        constexpr double roundFactor = 1 - std::numeric_limits<double>::epsilon();
        static_assert(roundFactor < 1);

        const int dim1 = int((factor * latt.row(0).norm() + roundFactor).getTrivial());
        const int dim2 = int((factor * latt.row(1).norm() + roundFactor).getTrivial());
        const int dim3 = int((factor * latt.row(2).norm() + roundFactor).getTrivial());
        return {dim1, dim2, dim3};
    }

    template<class ScalarType>
    ScalarType Ewald<ScalarType>::realSum(const LatticeMatrix& cellLattice,
                                          ScalarType integralLimit,
                                          std::tuple<int, int, int> dim,
                                          const Vector<ScalarType, 3>& deltaPos,
                                          const ScalarType& averageCellSize) {
        ScalarType sum = ScalarType::Zero();

        if (isTooClose(deltaPos.norm(), averageCellSize)) {
            for (int i = -std::get<0>(dim); i <= std::get<0>(dim); ++i) {
                for (int j = -std::get<1>(dim); j <= std::get<1>(dim); ++j) {
                    for (int k = -std::get<2>(dim); k <= std::get<2>(dim); ++k) {
                        if (i == 0 && j == 0 && k == 0)
                            continue;
                        const Vector3D transVector = Vector3D(cellLattice.row(0)) * ScalarType(i)
                                                   + Vector3D(cellLattice.row(1)) * ScalarType(j)
                                                   + Vector3D(cellLattice.row(2)) * ScalarType(k);
                        const ScalarType norm = transVector.norm();
                        sum += erfc(integralLimit * norm) / norm; //Optimize: VASP uses searching table method
                    }
                }
            }
        }
        else {
            for (int i = -std::get<0>(dim); i <= std::get<0>(dim); ++i) {
                for (int j = -std::get<1>(dim); j <= std::get<1>(dim); ++j) {
                    for (int k = -std::get<2>(dim); k <= std::get<2>(dim); ++k) {
                        const Vector3D transVector = Vector3D(cellLattice.row(0)) * ScalarType(i)
                                                   + Vector3D(cellLattice.row(1)) * ScalarType(j)
                                                   + Vector3D(cellLattice.row(2)) * ScalarType(k);
                        const ScalarType norm = (deltaPos - transVector).norm();
                        sum += erfc(integralLimit * norm) / norm;
                    }
                }
            }
        }
        return sum;
    }

    template<class ScalarType>
    ScalarType Ewald<ScalarType>::reciprocalSum(const LatticeMatrix& repLattice,
                                                ScalarType integralLimit,
                                                std::tuple<int, int, int> dim,
                                                const Vector<ScalarType, 3>& deltaPos) {
        ScalarType sum = ScalarType::Zero();
        for (int i = -std::get<0>(dim); i <= std::get<0>(dim); ++i) {
            for (int j = -std::get<1>(dim); j <= std::get<1>(dim); ++j) {
                for (int k = -std::get<2>(dim); k <= std::get<2>(dim); ++k) {
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

    template<class ScalarType>
    bool Ewald<ScalarType>::isTooClose(ScalarType deltaR, ScalarType averageCellSize) {
        return deltaR <= averageCellSize * std::numeric_limits<ScalarType>::epsilon();
    }
}
