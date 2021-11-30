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
#include "Physica/Core/Physics/ElectronStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronStructure/ReciprocalCell.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType>
        std::tuple<int, int, int> getSumDimention(const typename CrystalCell::LatticeMatrix& latt, ScalarType factor) {
            constexpr double roundFactor = 1 - std::numeric_limits<double>::epsilon();
            static_assert(roundFactor < 1);

            const int dim1 = int((factor * latt.row(0).norm() + roundFactor).getTrivial());
            const int dim2 = int((factor * latt.row(1).norm() + roundFactor).getTrivial());
            const int dim3 = int((factor * latt.row(2).norm() + roundFactor).getTrivial());
            return {dim1, dim2, dim3};
        }

        template<class ScalarType>
        ScalarType realSum(const CrystalCell& cell,
                           ScalarType integralLimit,
                           std::tuple<int, int, int> dim,
                           const Vector<ScalarType, 3>& deltaPos) {
            using VectorType = Vector<ScalarType, 3>;
            const auto& latt = cell.getLattice();
            ScalarType sum = ScalarType::Zero();

            if (deltaPos.squaredNorm().isZero()) {
                for (int i = -std::get<0>(dim); i <= std::get<0>(dim); ++i) {
                    for (int j = -std::get<1>(dim); j <= std::get<1>(dim); ++j) {
                        for (int k = -std::get<2>(dim); k <= std::get<2>(dim); ++k) {
                            if (i == 0 && j == 0 && k == 0)
                                continue;
                            const VectorType transVector = VectorType(latt.row(0)) * ScalarType(i)
                                                         + VectorType(latt.row(1)) * ScalarType(j)
                                                         + VectorType(latt.row(2)) * ScalarType(k);
                            const ScalarType norm = transVector.norm();
                            sum += erfc(integralLimit * norm) / norm;
                        }
                    }
                }
            }
            else {
                for (int i = -std::get<0>(dim); i <= std::get<0>(dim); ++i) {
                    for (int j = -std::get<1>(dim); j <= std::get<1>(dim); ++j) {
                        for (int k = -std::get<2>(dim); k <= std::get<2>(dim); ++k) {
                            const VectorType transVector = VectorType(latt.row(0)) * ScalarType(i)
                                                         + VectorType(latt.row(1)) * ScalarType(j)
                                                         + VectorType(latt.row(2)) * ScalarType(k);
                            const ScalarType norm = (deltaPos - transVector).norm();
                            sum += erfc(integralLimit * norm) / norm;
                        }
                    }
                }
            }
            return sum;
        }

        template<class ScalarType>
        ScalarType reciprocalSum(const ReciprocalCell& cell,
                                 ScalarType integralLimit,
                                 std::tuple<int, int, int> dim,
                                 const Vector<ScalarType, 3>& deltaPos) {
            ScalarType sum = ScalarType::Zero();
            for (int i = -std::get<0>(dim); i <= std::get<0>(dim); ++i) {
                for (int j = -std::get<1>(dim); j <= std::get<1>(dim); ++j) {
                    for (int k = -std::get<2>(dim); k <= std::get<2>(dim); ++k) {
                        if (i == 0 && j == 0 && k == 0)
                            continue;
                        const auto& latt = cell.getLattice();
                        using VectorType = Vector<ScalarType, 3>;
                        const VectorType repVector = VectorType(latt.row(0)) * ScalarType(i)
                                                   + VectorType(latt.row(1)) * ScalarType(j)
                                                   + VectorType(latt.row(2)) * ScalarType(k);
                        const ScalarType squaredNorm = repVector.squaredNorm();
                        const ScalarType dot = repVector * deltaPos;
                        sum += cos(dot) / (squaredNorm * exp(squaredNorm / square(ScalarType::Two() * integralLimit)));
                    }
                }
            }
            return sum;
        }
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     */
    template<class ScalarType>
    ScalarType getEwaldEnergy(const CrystalCell& cell, const ReciprocalCell& repCell) {
        const size_t ionCount = cell.getPos().getRow();
        const ScalarType inv_volume = reciprocal(cell.getVolume());
        const ScalarType averageCellSize = cbrt(ScalarType(cell.getVolume()));
        const auto realSumDim = Internal::getSumDimention(repCell.getLattice(), averageCellSize);
        const auto repSumDim = Internal::getSumDimention(cell.getLattice(), ScalarType(2 * M_PI) * reciprocal(averageCellSize));
        const ScalarType integralLimit = repCell.getMinNorm();

        ScalarType result = ScalarType::Zero();
        int totalCharge = 0;
        unsigned int totalSquaredCharge = 0;
        for (size_t ion1 = 0; ion1 < ionCount; ++ion1) {
            for (size_t ion2 = 0; ion2 < ionCount; ++ion2) {
                const Vector<ScalarType, 3> deltaPos = cell.getLattice() * (cell.getPos().row(ion1).asVector() - cell.getPos().row(ion2));
                ScalarType sum = Internal::realSum(cell, integralLimit, realSumDim, deltaPos);
                sum += ScalarType(4 * M_PI) * Internal::reciprocalSum(repCell, integralLimit, repSumDim, deltaPos) * inv_volume;

                const int charge1 = cell.getCharge(ion1);
                const int charge2 = cell.getCharge(ion2);
                const ScalarType dotCharge = ScalarType(charge1 * charge2);
                result += sum * dotCharge;
            }
            const int charge = cell.getCharge(ion1);
            totalCharge += charge;
            totalSquaredCharge += charge * charge;
        }
        result *= ScalarType(0.5);
        result -= ScalarType(totalSquaredCharge) * integralLimit / sqrt(ScalarType(M_PI));
        result -= ScalarType(totalCharge * totalCharge) * ScalarType(M_PI) / (ScalarType::Two() * square(integralLimit)) * inv_volume;
        return result;
    }
}
