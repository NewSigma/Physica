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
#include "PairModel.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType> class Q_TIP4P;

    namespace Internal {
        template<class T> class Traits;

        template<class ScalarType, class PosScalarType>
        class Traits<Q_TIP4P<ScalarType, PosScalarType>> {
        public:
            constexpr static double sigma = PhyConst<AU>::angstormToBohr(3.1589);

            static inline ScalarType lennardJonesPot(PosScalarType r) {
                constexpr double squared_sigma = sigma * sigma;
                const ScalarType rep_r2 = ScalarType(squared_sigma) / square(r);
                const ScalarType rep_r4 = square(rep_r2);
                const ScalarType rep_r6 = rep_r4 * rep_r2;
                const ScalarType rep_r12 = square(rep_r6);
                return rep_r12 - rep_r6;
            }

            static inline ScalarType lennardJonesForce(PosScalarType r) {
                const ScalarType rep_r = ScalarType(sigma) / r;
                const ScalarType rep_r2 = square(rep_r);
                const ScalarType rep_r4 = square(rep_r2);
                const ScalarType rep_r6 = rep_r4 * rep_r2;
                const ScalarType rep_r7 = rep_r6 * rep_r;
                const ScalarType rep_r13 = rep_r7 * rep_r6;
                return rep_r13 * 2 - rep_r7;
            }
        };
    }
    /**
     * Reference:
     * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
     * [2] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:205
     */
    template<class ScalarType, class PosScalarType>
    class Q_TIP4P final {
        using This = Q_TIP4P<ScalarType, PosScalarType>;
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using HytrogenListType = Utils::Array<std::pair<size_t, size_t>>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using EwaldType = Ewald<ScalarType, PosScalarType>;
        using LJModelType = PairModel<ScalarType, PosScalarType, decltype(&Internal::Traits<This>::lennardJonesPot)>;
        constexpr static unsigned int Dim = MDCellType::Dim;
    public:
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
        EwaldType ewald;
        ScalarType stepSize;
        LJModelType LJModel;
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
        [[nodiscard]] HytrogenListType& getHytrogenList() noexcept { return hytrogenList; }
        [[nodiscard]] size_t getNumMolecule() const noexcept { return hytrogenList.getLength(); }
        [[nodiscard]] const typename MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Helpers */
        void guessHytrogenList(const MDCellType& refer_cell);
        bool checkHytrogenList() const;
        void swap(Q_TIP4P& model) noexcept;
    private:
        Vector<ScalarType> makeCharges() const;
        PositionMatrix makeChargePos(const PositionMatrix& pos) const;
        ScalarType potentialEnergyWithoutEwald(const MDCellType& cell) const;
        ScalarType ewaldEnergy(const MDCellType& cell) const;
        ScalarType elasticEnergy(const MDCellType& cell) const;
        static ScalarType modifiedMorsePot(ScalarType r, size_t numMolecule);
        static ScalarType modifiedMorseForce(ScalarType r, size_t numMolecule);
    };

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>::Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_, ScalarType stepSize_)
            : hytrogenList(refer_cell.getNumParticle() / 3)
            , stepSize(std::move(stepSize_))
            , LJModel(std::move(cutoff_), Internal::Traits<This>::lennardJonesForce, Internal::Traits<This>::lennardJonesPot) {
        assert(refer_cell.getNumParticle() % 3 == 0);
        for (size_t i = 0; i < getNumMolecule(); ++i) {
            hytrogenList[i].first = i * 2;
            hytrogenList[i].second = i * 2 + 1;
        }
        ewald = EwaldType(refer_cell.getLattice(), makeCharges());
    }

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>& Q_TIP4P<ScalarType, PosScalarType>::operator=(Q_TIP4P model) noexcept {
        swap(model);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force(const MDCellType& cell) const {
        using Vector3D = Vector<PosScalarType, Dim>;
        const size_t numMolecule = getNumMolecule();
        Vector<ScalarType> result(3 * numMolecule * Dim, 0);
        /* LJ */ {
            const MDCellType cellWithoutH(cell.getLattice(), cell.getPos().bottomRows(2 * numMolecule), cell.getMassVec());
            const ScalarType factor = ScalarType((24 * epsilon / Internal::Traits<This>::sigma / PhyConst<SI>::avogadroNa) * numMolecule);
            auto force = result.tail(2 * numMolecule * Dim);
            force = LJModel.force(cellWithoutH) * factor;
        }
        /* Intra molecule */ {
            Vector3D vecOH1, vecOH2;
            Vector<ScalarType, 3> f;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                const size_t indexO = offset + i;
                const size_t indexH1 = hytrogenList[i].first;
                const size_t indexH2 = hytrogenList[i].second;
                vecOH1 = cell.minDistVector(indexO, indexH1);
                vecOH2 = cell.minDistVector(indexO, indexH2);
                const PosScalarType r1 = vecOH1.norm();
                const PosScalarType r2 = vecOH2.norm();
                auto forceO = result.segment(3 * indexO, 3 * indexO + 3);
                auto forceH1 = result.segment(3 * indexH1, 3 * indexH1 + 3);
                auto forceH2 = result.segment(3 * indexH2, 3 * indexH2 + 3);

                f = vecOH1 * (modifiedMorseForce(r1, numMolecule) / r1);
                forceO -= f;
                forceH1 += f;
                
                f = vecOH2 * (modifiedMorseForce(r2, numMolecule) / r2);
                forceO -= f;
                forceH1 += f;
            }

            for (size_t i = 0; i < result.getLength(); ++i) {
                result[i] += -Differential<ScalarType>::doublePoint([this, i, &cell](ScalarType x) -> ScalarType {
                    PositionMatrix pos = cell.getPos();
                    *(pos.begin() + i) = x;
                    MDCellType temp(cell.getLattice(), std::move(pos), cell.getMassVec());
                    return elasticEnergy(temp);
                }, ScalarType(cell.getPos().flatten().calc(i)), stepSize);
            }
        }
        Vector<ScalarType> coulomb;
        /* Coulomb Part */ {
            const PositionMatrix chargePos = makeChargePos(cell.getPos());
            coulomb = ewald.force(chargePos);

            PeriodicCell<PosScalarType, 3> chargeCell(cell.getLattice(), chargePos);
            const size_t minIndexO = 2 * numMolecule;
            const size_t maxIndexO = minIndexO + numMolecule;
            Vector<ScalarType, 3> f;
            Vector3D vecMH1, vecMH2, vecH1H2;
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                const size_t indexH1 = hytrogenList[i - minIndexO].first;
                const size_t indexH2 = hytrogenList[i - minIndexO].second;
                auto forceM = coulomb.segment(3 * i, 3 * i + 3);
                auto forceH1 = coulomb.segment(3 * indexH1, 3 * indexH1 + 3);
                auto forceH2 = coulomb.segment(3 * indexH2, 3 * indexH2 + 3);

                vecMH1 = chargeCell.minDistVector(i, indexH1);
                PosScalarType r2 = vecMH1.squaredNorm();
                f = vecMH1 * (ScalarType(-charge * charge * 0.5) / (r2 * sqrt(r2)));
                forceM += f;
                forceH1 -= f;

                vecMH2 = chargeCell.minDistVector(i, indexH2);
                r2 = vecMH2.squaredNorm();
                f = vecMH2 * (ScalarType(-charge * charge * 0.5) / (r2 * sqrt(r2)));
                forceM += f;
                forceH2 -= f;

                vecH1H2 = chargeCell.minDistVector(indexH1, indexH2);
                r2 = vecH1H2.squaredNorm();
                f = vecH1H2 * (ScalarType(charge * charge * 0.25) / (r2 * sqrt(r2)));
                forceH1 += f;
                forceH2 -= f;
            }
            /* Change representation: from partial charge representation to HOH representation */
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                const size_t indexH1 = hytrogenList[i - minIndexO].first;
                const size_t indexH2 = hytrogenList[i - minIndexO].second;
                auto forceO = coulomb.segment(3 * i, 3 * i + 3);
                auto forceH1 = coulomb.segment(3 * indexH1, 3 * indexH1 + 3);
                auto forceH2 = coulomb.segment(3 * indexH2, 3 * indexH2 + 3);
                f = forceO * ScalarType((1 - gamma) * 0.5);
                forceH1 += f;
                forceH2 += f;
                forceO *= ScalarType(gamma);
            }
        }
        return result + coulomb;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergy(const MDCellType& cell) const {
        return potentialEnergyWithoutEwald(cell) + ewaldEnergy(cell);
    }

    template<class ScalarType, class PosScalarType>
    void Q_TIP4P<ScalarType, PosScalarType>::guessHytrogenList(const MDCellType& refer_cell) {
        const auto& pos = refer_cell.getPos();
        const size_t numMolecule = getNumMolecule();
        const size_t maxHytrogenId = refer_cell.getNumParticle() - numMolecule;
        ScalarType bondLength1, bondLength2;
        size_t hytrogenId1 = 0, hytrogenId2 = 0;
        for (size_t i = 0; i < numMolecule; ++i) {
            bondLength1 = bondLength2 = std::numeric_limits<ScalarType>::max();

            auto pos_O = pos.row(maxHytrogenId + i);
            for (size_t j = 0; j < maxHytrogenId; ++j) {
                auto pos_H = pos.row(j);
                const ScalarType squared_dist = refer_cell.minDistVector(maxHytrogenId + i, j).squaredNorm();
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
        ewald.swap(model.ewald);
        stepSize.swap(model.stepSize);
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::makeCharges() const {
        const size_t numMolecule = getNumMolecule();
        const size_t maxIndexH = 2 * numMolecule;
        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType> charges(maxIndexO);
        for (size_t i = 0; i < maxIndexH; ++i)
            charges[i] = ScalarType(charge * 0.5);
        
        for (size_t i = maxIndexH; i < maxIndexO; ++i)
            charges[i] = ScalarType(-charge);
        return charges;
    }

    template<class ScalarType, class PosScalarType>
    typename Q_TIP4P<ScalarType, PosScalarType>::PositionMatrix
    Q_TIP4P<ScalarType, PosScalarType>::makeChargePos(const PositionMatrix& pos) const {
        PositionMatrix chargePos(pos.getRow(), 3);
        const size_t numMolecule = pos.getRow() / 3;
        const size_t maxIndexH = 2 * numMolecule;
        auto chargePosH = chargePos.topRows(maxIndexH);
        chargePosH = pos.topRows(maxIndexH);

        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numMolecule;
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            auto chargePosO = chargePos.row(i);
            auto posO = pos.row(i);
            auto posH1 = pos.row(hytrogenList[i - minIndexO].first);
            auto posH2 = pos.row(hytrogenList[i - minIndexO].second);
            chargePosO = posO.asVector() * ScalarType(gamma) + (posH1.asVector() + posH2.asVector()) * ScalarType((1 - gamma) * 0.5);
        }
        return chargePos;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergyWithoutEwald(const MDCellType& cell) const {
        using Vector3D = Vector<PosScalarType, Dim>;
        const size_t numMolecule = getNumMolecule();

        const MDCellType cellWithoutH(cell.getLattice(), cell.getPos().bottomRows(2 * numMolecule), cell.getMassVec());
        const ScalarType factor = (4 * epsilon / PhyConst<SI>::avogadroNa) * getNumMolecule();
        const ScalarType interMoleculeEnergy = LJModel.potentialEnergy(cellWithoutH) * factor;
        ScalarType intraMoleculeEnergy = 0;
        /* Intra molecule */ {
            Vector3D vecOH1, vecOH2;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                vecOH1 = cell.minDistVector(offset + i, hytrogenList[i].first);
                vecOH2 = cell.minDistVector(offset + i, hytrogenList[i].second);
                const ScalarType r1 = vecOH1.norm();
                const ScalarType r2 = vecOH2.norm();
                const ScalarType elastic = square(arccos((vecOH1 * vecOH2) / (r1 * r2)) - ScalarType(equalTheta)) * (kTheta / PhyConst<SI>::avogadroNa * numMolecule * 0.5);
                intraMoleculeEnergy += modifiedMorsePot(r1, numMolecule) + modifiedMorsePot(r2, numMolecule) + elastic;
            }
        }

        return interMoleculeEnergy + intraMoleculeEnergy;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::ewaldEnergy(const MDCellType& cell) const {
        const PositionMatrix chargePos = makeChargePos(cell.getPos());
        const size_t numMolecule = getNumMolecule();
        ScalarType selfCoulomb = 0;
        /* Spurious Coulomb Part */ {
            PeriodicCell<PosScalarType, 3> chargeCell(cell.getLattice(), chargePos);
            const size_t minIndexO = 2 * numMolecule;
            const size_t maxIndexO = minIndexO + numMolecule;
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                const size_t indexH1 = hytrogenList[i - minIndexO].first;
                const size_t indexH2 = hytrogenList[i - minIndexO].second;
                selfCoulomb += ScalarType(-charge * charge * 0.5) / chargeCell.minDistVector(i, indexH1).norm();
                selfCoulomb += ScalarType(-charge * charge * 0.5) / chargeCell.minDistVector(i, indexH2).norm();
                selfCoulomb += ScalarType(charge * charge * 0.25) / chargeCell.minDistVector(indexH1, indexH2).norm();
            }
        }
        return ewald.potentialEnergy(chargePos) - selfCoulomb;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::elasticEnergy(const MDCellType& cell) const {
        using Vector3D = Vector<PosScalarType, Dim>;
        const size_t numMolecule = getNumMolecule();
        ScalarType result = 0;
        Vector3D vecOH1, vecOH2;
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            vecOH1 = cell.minDistVector(offset + i, hytrogenList[i].first);
            vecOH2 = cell.minDistVector(offset + i, hytrogenList[i].second);
            const ScalarType r1 = vecOH1.norm();
            const ScalarType r2 = vecOH2.norm();
            result += square(arccos((vecOH1 * vecOH2) / (r1 * r2)) - ScalarType(equalTheta)) * (kTheta / PhyConst<SI>::avogadroNa * numMolecule * 0.5);
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::modifiedMorsePot(ScalarType r, size_t numMolecule) {
        const ScalarType delta = (r - ScalarType(equalR)) * alphaR;
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        const ScalarType delta4 = square(delta2);
        return (delta2 - delta3 + delta4 * (7.0 / 12)) * ((Dr / PhyConst<SI>::avogadroNa) * numMolecule);
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::modifiedMorseForce(ScalarType r, size_t numMolecule) {
        const ScalarType delta = (r - ScalarType(equalR)) * alphaR;
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        const ScalarType delta4 = square(delta2);
        return -(delta2 * 2 - delta3 * 3 + delta4 * (7.0 / 3)) * ((Dr * alphaR / PhyConst<SI>::avogadroNa) * numMolecule);
    }
}
