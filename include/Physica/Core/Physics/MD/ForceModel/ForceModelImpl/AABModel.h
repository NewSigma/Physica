/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermutationMatrix.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    /**
     * AAB(2 : 1) molecular, for example H2O(\class Q_TIP4P)
     */
    template<class ScalarType>
    class AABModel {
    public:
        constexpr static unsigned int Dim = 3;
        using PlainScalar = typename ScalarType::PlainScalar;
        using MDCellType = MDCell<ScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
    public:
        ~AABModel() = default;
        /* Static members */
        static PermutationMatrix<ScalarType> sortPosition(MDCellType& cell, size_t atomicNum1, size_t atomicNum2);
        static bool isCellOrdered(const MDCellType& cell, size_t atomicNum1, size_t atomicNum2);
        static Vector<ScalarType> makeCharges(size_t numMolecule, PlainScalar atomCharge1, PlainScalar atomCharge2);
    protected:
        AABModel() = default;
        AABModel(const AABModel&) = default;
        AABModel(AABModel&&) noexcept = default;
        /* Operators */
        AABModel& operator=(const AABModel&) = default;
        AABModel& operator=(AABModel&&) noexcept = default;
    };

    template<class ScalarType>
    PermutationMatrix<ScalarType> AABModel<ScalarType>::sortPosition(MDCellType& cell, size_t atomicNum1, size_t atomicNum2) {
        using MassVector = typename MDCellType::MassVector;
        const auto& source = cell.getPos();
        const size_t numAtom = source.getRow();
        assert(numAtom % 3 == 0);
        const size_t numH = numAtom * 2 / 3;
        const size_t numO = numAtom / 3;

        PositionMatrix new_pos(source.getRow(), 3);
        MassVector new_mass(numAtom);
        Utils::Array<size_t> orderStage1(numAtom);
        /* Stage 1: Classify H and O */ {
            size_t indexH = 0, indexO = numH;
            for (size_t i = 0; i < numAtom; ++i) {
                const bool isHydrogen = cell.getMass(i).getValue() == PhyConst<AU>::atomMass(atomicNum1);
                const size_t index = isHydrogen ? indexH : indexO;
                new_pos.row(index) = source.row(i);
                new_mass[i] = i < numH ? PhyConst<AU>::atomMass(atomicNum1) : PhyConst<AU>::atomMass(atomicNum2);
                orderStage1[index] = i;
                indexH += isHydrogen;
                indexO += !isHydrogen;
            }
            assert(indexH == numH);
            assert(indexO == numAtom);
            cell.setPos(new_pos);
        }
        Utils::Array<size_t> orderStage2(numAtom);
        /* Stage 2: Sort H */ {
            for (size_t i = 0; i < numO; ++i) {
                const size_t indexO = i + numH;
                size_t indexH1 = 0, indexH2 = 0;
                /* Make indexH1, indexH2 */ {
                    ScalarType dist1, dist2;
                    dist1 = dist2 = std::numeric_limits<ScalarType>::max();
                    
                    for (size_t j = 0; j < numH; ++j) {
                        auto posOH = cell.minDistVector(indexO, j);
                        const ScalarType dist = posOH.squaredNorm();
                        if (dist1 > dist2) {
                            if (dist1 > dist) {
                                dist1 = dist;
                                indexH1 = j;
                            }
                        }
                        else {
                            if (dist2 > dist) {
                                dist2 = dist;
                                indexH2 = j;
                            }
                        }
                    }
                    if (indexH1 > indexH2)
                        std::swap(indexH1, indexH2);
                }
                auto posH1 = new_pos.row(2 * i);
                posH1 = source.row(indexH1).asVector();
                orderStage2[2 * i] = orderStage1[indexH1];

                auto posH2 = new_pos.row(2 * i + 1);
                posH2 = source.row(indexH2).asVector();
                orderStage2[2 * i + 1] = orderStage1[indexH2];
                orderStage2[indexO] = orderStage1[indexO];
            }
            cell = MDCellType(cell.getLattice(), std::move(new_pos), std::move(new_mass));
        }
        return PermutationMatrix<ScalarType>(std::move(orderStage2));
    }

    template<class ScalarType>
    bool AABModel<ScalarType>::isCellOrdered(const MDCellType& cell, size_t atomicNum1, size_t atomicNum2) {
        const size_t numParticle= cell.getNumParticle();
        const size_t maxIndexH = 2 * numParticle / 3;
        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numParticle / 3;
        for (size_t i = 0; i < maxIndexH; ++i) {
            const bool isHydrogen = cell.getMass(i).getValue() == PhyConst<AU>::atomMass(atomicNum1);
            if (!isHydrogen)
                return false;
        }
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            const bool isOxygen = cell.getMass(i).getValue() == PhyConst<AU>::atomMass(atomicNum2);
            if (!isOxygen)
                return false;
        }
        return true;
    }

    template<class ScalarType>
    Vector<ScalarType> AABModel<ScalarType>::makeCharges(size_t numMolecule, PlainScalar atomCharge1, PlainScalar atomCharge2) {
        const size_t maxIndexH = 2 * numMolecule;
        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType> charges(maxIndexO);
        for (size_t i = 0; i < maxIndexH; ++i)
            charges[i] = atomCharge1;
        
        for (size_t i = minIndexO; i < maxIndexO; ++i)
            charges[i] = atomCharge2;
        return charges;
    }
}
