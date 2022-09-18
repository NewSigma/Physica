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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Flatten.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/Ewald.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
     */
    template<class ScalarType, class PosScalarType>
    class Q_TIP4P {
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        constexpr static unsigned int Dim = MDCellType::Dim;
        using HytrogenListType = Utils::Array<std::pair<size_t, size_t>>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using EwaldType = Ewald<ScalarType, PosScalarType>;

        constexpr static double sigma = PhyConst<AU>::angstormToBohr(3.1589);
        constexpr static double epsilon = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(0.1852 * 1000) / PhyConst<SI>::unitCharge);
        constexpr static double charge = 1.1128;
        constexpr static double gamma = 0.73612;
        constexpr static double Dr = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(116.9 * 1000) / PhyConst<SI>::unitCharge);
        constexpr static double alphaR = PhyConst<AU>::bohrToAngstorm(2.287);
        constexpr static double equalR = PhyConst<AU>::angstormToBohr(0.9419);
        constexpr static double kTheta = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(87.85 * 1000) / PhyConst<SI>::unitCharge);
        constexpr static double equalTheta = PhyConst<SI>::degreeToRadian(107.4);
    private:
        HytrogenListType hytrogenList;
        PositionMatrix chargePos;
        EwaldType ewald;
        ScalarType cutoff;
        ScalarType stepSize;
    public:
        Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_, ScalarType stepSize_);
        Q_TIP4P(const Q_TIP4P&) = default;
        Q_TIP4P(Q_TIP4P&&) noexcept = default;
        ~Q_TIP4P() = default;
        /* Operators */
        Q_TIP4P& operator=(Q_TIP4P model) noexcept;
        /* Operations */
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        /* Getters */
        [[nodiscard]] const HytrogenListType& getHytrogenList() const noexcept { return hytrogenList; }
        [[nodiscard]] size_t getNumMolecule() const noexcept { return hytrogenList.getLength(); }
        /* Helpers */
        bool checkHytrogenList() const;
        void swap(Q_TIP4P& model) noexcept;
    private:
        void makeHytrogenList(const MDCellType& refer_cell);
        Vector<ScalarType> makeCharges() const;
        ScalarType minSquaredDist(const MDCellType& cell, size_t from_id, size_t to_id);
        ScalarType modifiedMorse(ScalarType r, size_t numMolecule) const;
        ScalarType potentialEnergyWithoutCoulomb(const MDCellType& cell) const;
    };

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>::Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_, ScalarType stepSize_)
            : chargePos(refer_cell.getNumParticle(), 3)
            , cutoff(std::move(cutoff_))
            , stepSize(std::move(stepSize_)) {
        assert(refer_cell.getNumParticle() % 3 == 0);
        makeHytrogenList(refer_cell);
        ewald = EwaldType(refer_cell.getLattice(), makeCharges());
    }

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>& Q_TIP4P<ScalarType, PosScalarType>::operator=(Q_TIP4P model) noexcept {
        swap(model);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force(const MDCellType& cell) const {
        Vector<ScalarType> result(getNumMolecule() * Dim);
        for (size_t i = 0; i < result.getLength(); ++i) {
            result[i] = -Differential<ScalarType>::doublePoint([this, i, &cell](ScalarType x) -> ScalarType {
                PositionMatrix pos = cell.getPos();
                *(pos.begin() + i) = x;
                return potentialEnergyWithoutCoulomb(MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec()));
            }, ScalarType(cell.getPos().flatten().calc(i)), stepSize);
        }
        return result + ewald.force(cell.getPos());
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergy(const MDCellType& cell) const {
        const size_t numMolecule = getNumMolecule();
        const auto& pos = cell.getPos();

        ScalarType coulomb = 0;
        /* Coulomb Part */ {
            const size_t maxIndexH = 2 * numMolecule;
            auto posH = chargePos.topRows(maxIndexH);
            posH = pos.topRows(maxIndexH);

            ScalarType selfCoulomb = 0;
            const size_t minIndexO = maxIndexH;
            const size_t maxIndexO = minIndexO + numMolecule;
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                auto partialCharge = chargePos.row(i);
                auto posO = pos.row(i);
                auto posH1 = pos.row(hytrogenList[i - minIndexO].first);
                auto posH2 = pos.row(hytrogenList[i - minIndexO].second);
                partialCharge = posO.asVector() * ScalarType(gamma) + (posH1.asVector() + posH2) * ScalarType((1 - gamma) * 0.5);

                selfCoulomb += ScalarType(-charge * charge / 2) / ScalarType((posH1.asVector() - partialCharge).norm());
                selfCoulomb += ScalarType(-charge * charge / 2) / ScalarType((posH2.asVector() - partialCharge).norm());
            }
            coulomb = ewald(pos) - selfCoulomb;
        }
        return potentialEnergyWithoutCoulomb(cell) + coulomb;
    }

    template<class ScalarType, class PosScalarType>
    bool Q_TIP4P<ScalarType, PosScalarType>::checkHytrogenList() const {
        const size_t maxHytrogenId = hytrogenList.getLength() * 2;
        for (size_t i = 0; i < hytrogenList.getLength(); ++i) {
            const auto& pair = hytrogenList[i];
            if (pair.first == pair.second)
                return false;

            if (pair.first >= maxHytrogenId || pair.second >= maxHytrogenId)
                return false;

            for (size_t j = i + 1; j < hytrogenList.getLength(); ++j) {
                const auto& pair2 = hytrogenList[j];
                if (pair.first == pair2.first
                 || pair.first == pair2.second
                 || pair.second == pair2.first
                 || pair.second == pair2.second)
                    return false;
            }
        }
        return true;
    }

    template<class ScalarType, class PosScalarType>
    void Q_TIP4P<ScalarType, PosScalarType>::swap(Q_TIP4P& model) noexcept {
        hytrogenList.swap(model.hytrogenList);
        chargePos.swap(model.chargePos);
        ewald.swap(model.ewald);
        cutoff.swap(model.cutoff);
        stepSize.swap(model.stepSize);
    }

    template<class ScalarType, class PosScalarType>
    void Q_TIP4P<ScalarType, PosScalarType>::makeHytrogenList(const MDCellType& refer_cell) {
        const size_t numMolecule = refer_cell.getNumParticle() / 3;
        hytrogenList.resize(numMolecule);

        const auto& pos = refer_cell.getPos();
        const size_t maxHytrogenId = refer_cell.getNumParticle() - numMolecule;
        ScalarType bondLength1, bondLength2;
        size_t hytrogenId1 = 0, hytrogenId2 = 0;
        for (size_t i = 0; i < numMolecule; ++i) {
            bondLength1 = bondLength2 = std::numeric_limits<ScalarType>::max();

            auto pos_O = pos.row(maxHytrogenId + i);
            for (size_t j = 0; j < maxHytrogenId; ++j) {
                auto pos_H = pos.row(j);
                const ScalarType squared_dist = minSquaredDist(refer_cell, maxHytrogenId + i, j);
                if (squared_dist < bondLength1) {
                    if (bondLength2 > bondLength1) {
                        bondLength2 = squared_dist;
                        hytrogenId2 = j;
                    }
                    else {
                        bondLength1 = squared_dist;
                        hytrogenId1 = j;
                    }
                }
                else if (squared_dist < bondLength2) {
                    bondLength2 = squared_dist;
                    hytrogenId2 = j;
                }
            }
            auto& hytrogenPair = hytrogenList[i];
            hytrogenPair.first = hytrogenId1;
            hytrogenPair.second = hytrogenId2;
        }
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::makeCharges() const {
        const size_t numMolecule = getNumMolecule();
        const size_t maxIndexH = 2 * numMolecule;
        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType> charges(maxIndexO);
        for (size_t i = 0; i < maxIndexH; ++i)
            charges[i] = ScalarType(-charge * 0.5);
        
        for (size_t i = maxIndexH; i < maxIndexO; ++i)
            charges[i] = ScalarType(charge);
        return charges;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::minSquaredDist(const MDCellType& cell, size_t from_id, size_t to_id) {
        using Vector3D = Vector<ScalarType, 3>;

        const auto& pos = cell.getPos();
        const auto& lattice = cell.getLattice();
        auto pos1 = pos.row(from_id);
        auto pos2 = pos.row(to_id);

        ScalarType result = std::numeric_limits<ScalarType>::max();
        for (int x = -1; x <= 1; ++x) {
            const Vector3D v1 = pos2.asVector() + lattice.row(0).asVector() * ScalarType(x);
            for (int y = -1; y <= 1; ++y) {
                const Vector3D v2 = v1 + lattice.row(1).asVector() * ScalarType(y);
                for (int z = -1; z <= 1; ++z) {
                    const Vector3D v3 = v2 + lattice.row(2).asVector() * ScalarType(z);
                    result = std::min(result, (v3 - pos1.asVector()).squaredNorm());
                }
            }
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::modifiedMorse(ScalarType r, size_t numMolecule) const {
        const ScalarType delta = (r - ScalarType(equalR)) * alphaR;
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        const ScalarType delta4 = square(delta2);
        return (delta2 - delta3 + delta4 * (7.0 / 12)) * ((Dr / PhyConst<SI>::avogadroNa) * numMolecule);
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergyWithoutCoulomb(const MDCellType& cell) const {
        using VectorType = Vector<ScalarType, Dim>;

        const size_t numMolecule = getNumMolecule();
        const auto& lattice = cell.getLattice();
        const auto& pos = cell.getPos();
        const auto range = MDCellType::estimateRange(lattice, cutoff);
        const double epsilonE4 = (4 * epsilon / PhyConst<SI>::avogadroNa) * numMolecule;

        ScalarType interMoleculeEnergy = 0;
        /* Inter molecular */ {
            /* LJ Part */
            MDCellType::forParticleInRange(range, lattice,
                [this, numMolecule, pos, epsilonE4, &interMoleculeEnergy](VectorType delta) {
                    ScalarType energy = 0;
                    VectorType posO1, posH11, posH12;
                    const size_t offset = 2 * numMolecule;
                    for (size_t i = 0; i < numMolecule; ++i) {
                        posO1 = pos.row(offset + i) + delta;
                        for (size_t j = i; j < numMolecule; ++j) {
                            auto posO2 = pos.row(offset + j);
                            const ScalarType r2 = (posO1 - posO2).squaredNorm();
                            const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                            if (isNotSelf && r2 < square(cutoff)) {
                                constexpr double squared_sigma = sigma * sigma;
                                const ScalarType rep_r2 = reciprocal(r2) * ScalarType(squared_sigma);
                                const ScalarType rep_r4 = square(rep_r2);
                                const ScalarType rep_r6 = rep_r4 * rep_r2;
                                const ScalarType rep_r12 = square(rep_r6);
                                energy += (rep_r12 - rep_r6) * epsilonE4;
                            }
                        }
                    }
                    interMoleculeEnergy += energy;
                });
        }
        ScalarType intraMoleculeEnergy = 0;
        /* Intra molecular */ {
            VectorType vecOH1, vecOH2;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                auto posO = pos.row(offset + i);
                auto posH1 = pos.row(hytrogenList[i].first);
                auto posH2 = pos.row(hytrogenList[i].second);

                vecOH1 = posH1.asVector() - posO;
                vecOH2 = posH2.asVector() - posO;
                const ScalarType r1 = vecOH1.norm();
                const ScalarType r2 = vecOH2.norm();
                const ScalarType elastic = square(arccos((vecOH1 * vecOH2) / (r1 * r2)) - ScalarType(equalTheta)) * (kTheta * 0.5);
                intraMoleculeEnergy += modifiedMorse(r1, numMolecule) + modifiedMorse(r2, numMolecule);
            }
        }
        return interMoleculeEnergy + intraMoleculeEnergy;
    }
}
