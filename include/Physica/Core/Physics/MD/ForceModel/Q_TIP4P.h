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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermutationMatrix.h"
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
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using EwaldType = Ewald<ScalarType, PosScalarType>;
        using LJModelType = PairModel<ScalarType, PosScalarType, decltype(&Internal::Traits<This>::lennardJonesPot)>;
        constexpr static unsigned int Dim = MDCellType::Dim;
    public:
        constexpr static double epsilon = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(0.1852 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double charge = 1.1128;
        constexpr static double gamma = 0.73612;
        constexpr static double Dr = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(116.09 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double alphaR = PhyConst<AU>::bohrToAngstorm(2.287);
        constexpr static double equalR = PhyConst<AU>::angstormToBohr(0.9419);
        constexpr static double kTheta = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(87.85 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double equalTheta = PhyConst<SI>::degreeToRadian(107.4);
    private:
        size_t numMolecule;
        EwaldType ewald;
        LJModelType LJModel;
    public:
        Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_);
        Q_TIP4P(const Q_TIP4P&) = default;
        Q_TIP4P(Q_TIP4P&&) noexcept = default;
        ~Q_TIP4P() = default;
        /* Operators */
        Q_TIP4P& operator=(Q_TIP4P model) noexcept;
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const;
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return numMolecule; }
        [[nodiscard]] const typename MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Helpers */
        void swap(Q_TIP4P& model) noexcept;
        /* Static members */
        [[nodiscard]] static PositionMatrix makeDipoleMoments(const PeriodicCell<PosScalarType, 3>& cell);
        static PermutationMatrix<PosScalarType> sortPosition(PeriodicCell<PosScalarType, 3>& cell);
    private:
        Vector<ScalarType> makeCharges() const;
        PositionMatrix makeChargePos(const MDCellType& cell) const;
        ScalarType potentialEnergyWithoutEwald(const MDCellType& cell) const;
        ScalarType ewaldEnergy(const MDCellType& cell) const;
        static ScalarType modifiedMorsePot(ScalarType r);
        static ScalarType modifiedMorseForce(ScalarType r);
    };

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>::Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_)
            : numMolecule(refer_cell.getNumParticle() / 3)
            , LJModel(std::move(cutoff_), Internal::Traits<This>::lennardJonesForce, Internal::Traits<This>::lennardJonesPot) {
        assert(refer_cell.getNumParticle() % 3 == 0);
        ewald = EwaldType(refer_cell.getLattice(), makeCharges());
    }

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>& Q_TIP4P<ScalarType, PosScalarType>::operator=(Q_TIP4P model) noexcept {
        swap(model);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0);

        Vector<ScalarType> result;
        auto future = Executor::schedule([this, &cell, &result]() {
            result = force_short<Executor>(cell);
        });

        const Vector<ScalarType> coulomb = force_long<Executor>(cell);
        Executor::auto_wait(future);
        result += coulomb;
        return result;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force_short(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0);
        using Vector3D = Vector<PosScalarType, Dim>;

        Vector<ScalarType> shortForce(3 * numMolecule * Dim, 0);
        /* LJ */ {
            const MDCellType cellWithoutH(cell.getLattice(), cell.getPos().bottomRows(2 * numMolecule), cell.getMassVec());
            const ScalarType factor = ScalarType(24 * epsilon / Internal::Traits<This>::sigma);
            auto force = shortForce.tail(2 * numMolecule * Dim);
            force = LJModel.template force<Executor>(cellWithoutH) * factor;
        }
        /* Intra molecule */ {
            Vector3D vecOH1, vecOH2;
            Vector<ScalarType, 3> f;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                const size_t indexO = offset + i;
                const size_t indexH1 = 2 * i;
                const size_t indexH2 = 2 * i + 1;
                vecOH1 = cell.minDistVector(indexO, indexH1);
                vecOH2 = cell.minDistVector(indexO, indexH2);
                const PosScalarType r1 = vecOH1.norm();
                const PosScalarType r2 = vecOH2.norm();
                auto forceO = shortForce.segment(3 * indexO, 3 * indexO + 3);
                auto forceH1 = shortForce.segment(3 * indexH1, 3 * indexH1 + 3);
                auto forceH2 = shortForce.segment(3 * indexH2, 3 * indexH2 + 3);

                f = vecOH1 * (modifiedMorseForce(r1) / r1);
                forceO -= f;
                forceH1 += f;
                    
                f = vecOH2 * (modifiedMorseForce(r2) / r2);
                forceO -= f;
                forceH2 += f;
            }

            Vector3D force1, force2;
            for (size_t i = 0; i < numMolecule; ++i) {
                vecOH1 = cell.minDistVector(offset + i, 2 * i);
                vecOH2 = cell.minDistVector(offset + i, 2 * i + 1);
                const PosScalarType rep_r1 = reciprocal(vecOH1.norm());
                const PosScalarType rep_r2 = reciprocal(vecOH2.norm());
                const ScalarType u = (vecOH1 * vecOH2) * (rep_r1 * rep_r2);
                const ScalarType temp = reciprocal(sqrt(ScalarType(1) - square(u)));
                const ScalarType factor1 = square(rep_r1) * u * temp;
                const ScalarType factor2 = square(rep_r2) * u * temp;
                const ScalarType factor3 = rep_r1 * rep_r2 * temp;

                force1 = vecOH1 * factor1 - vecOH2 * factor3;
                force2 = vecOH2 * factor2 - vecOH1 * factor3;
                const ScalarType theta = arccos(u);
                const ScalarType factor = (ScalarType(equalTheta) - theta) * kTheta;
                force1 *= factor;
                force2 *= factor;

                const size_t indexO = offset + i;
                const size_t indexH1 = 2 * i;
                const size_t indexH2 = 2 * i + 1;
                auto forceO = shortForce.segment(3 * indexO, 3 * indexO + 3);
                auto forceH1 = shortForce.segment(3 * indexH1, 3 * indexH1 + 3);
                auto forceH2 = shortForce.segment(3 * indexH2, 3 * indexH2 + 3);
                forceH1 += force1;
                forceH2 += force2;
                forceO -= (force1 + force2);
            }
        }
        return shortForce;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force_long(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0);
        using Vector3D = Vector<PosScalarType, Dim>;

        Vector<ScalarType> coulomb;
        /* Coulomb Part */ {
            const PositionMatrix chargePos = makeChargePos(cell);
            coulomb = ewald.template force<Executor>(chargePos);
            PeriodicCell<PosScalarType, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
            const size_t minIndexO = 2 * numMolecule;
            const size_t maxIndexO = minIndexO + numMolecule;
            Vector<ScalarType, 3> f;
            Vector3D vecMH1, vecMH2, vecH1H2;
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                const size_t indexH1 = 2 * (i - minIndexO);
                const size_t indexH2 = indexH1 + 1;
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
                const size_t indexH1 = 2 * (i - minIndexO);
                const size_t indexH2 = indexH1 + 1;
                auto forceO = coulomb.segment(3 * i, 3 * i + 3);
                auto forceH1 = coulomb.segment(3 * indexH1, 3 * indexH1 + 3);
                auto forceH2 = coulomb.segment(3 * indexH2, 3 * indexH2 + 3);
                f = forceO * ScalarType((1 - gamma) * 0.5);
                forceH1 += f;
                forceH2 += f;
                forceO *= ScalarType(gamma);
            }
        }
        return coulomb;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergy(const MDCellType& cell) const {
        return potentialEnergyWithoutEwald(cell) + ewaldEnergy(cell);
    }

    template<class ScalarType, class PosScalarType>
    void Q_TIP4P<ScalarType, PosScalarType>::swap(Q_TIP4P& model) noexcept {
        std::swap(numMolecule, model.numMolecule);
        ewald.swap(model.ewald);
    }

    template<class ScalarType, class PosScalarType>
    typename Q_TIP4P<ScalarType, PosScalarType>::PositionMatrix
    Q_TIP4P<ScalarType, PosScalarType>::makeDipoleMoments(const PeriodicCell<PosScalarType, 3>& cell) {
        assert(cell.getNumParticle() % 3 == 0);
        const size_t numMolecule = cell.getNumParticle() / 3;
        PositionMatrix dipoles(numMolecule, 3);
        for (size_t i = 0; i < dipoles.getRow(); ++i) {
            const size_t indexO = 2 * numMolecule + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            auto dipole = dipoles.row(i);
            dipole = (cell.minDistVector(indexO, indexH1) + cell.minDistVector(indexO, indexH2)) * ScalarType(-gamma * charge);
        }
        return dipoles;
    }

    template<class ScalarType, class PosScalarType>
    PermutationMatrix<PosScalarType> Q_TIP4P<ScalarType, PosScalarType>::sortPosition(PeriodicCell<PosScalarType, 3>& cell) {
        const auto& source = cell.getPos();
        const size_t numAtom = source.getRow();
        assert(numAtom % 3 == 0);
        const size_t numH = numAtom * 2 / 3;
        const size_t numO = numAtom / 3;

        PositionMatrix new_pos(source.getRow(), 3);
        Utils::Array<size_t> order(numAtom);
        for (size_t i = 0; i < numO; ++i) {
            const size_t indexO = i + numH;
            auto posO = new_pos.row(indexO);
            posO = source.row(indexO).asVector();

            size_t indexH1 = 0, indexH2 = 0;
            /* Make indexH1, indexH2 */ {
                PosScalarType dist1, dist2;
                dist1 = dist2 = std::numeric_limits<PosScalarType>::max();
                
                for (size_t j = 0; j < numH; ++j) {
                    auto posOH = cell.minDistVector(indexO, j);
                    const PosScalarType dist = posOH.squaredNorm();
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
            order[2 * i] = indexH1;

            auto posH2 = new_pos.row(2 * i + 1);
            posH2 = source.row(indexH2).asVector();
            order[2 * i + 1] = indexH2;
            order[indexO] = indexO;
        }
        cell.setPos(new_pos);
        return PermutationMatrix<PosScalarType>(std::move(order));
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
    Q_TIP4P<ScalarType, PosScalarType>::makeChargePos(const MDCellType& cell) const {
        PositionMatrix chargePos(cell.getPos().getRow(), 3);
        const auto& pos = cell.getPos();
        const size_t numMolecule = getNumMolecule();
        const size_t maxIndexH = 2 * numMolecule;
        auto chargePosH = chargePos.topRows(maxIndexH);
        chargePosH = pos.topRows(maxIndexH);

        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<PosScalarType, 3> vecOH1, vecOH2;
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            auto chargePosO = chargePos.row(i);
            auto posO = pos.row(i);
            const size_t indexH1 = 2 * (i - minIndexO);
            const size_t indexH2 = indexH1 + 1;
            vecOH1 = cell.minDistVector(i, indexH1);
            vecOH2 = cell.minDistVector(i, indexH2);
            chargePosO = posO.asVector() + (vecOH1 + vecOH2) * ScalarType((1 - gamma) * 0.5);
        }
        cell.normalizePos(chargePos);
        return chargePos;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergyWithoutEwald(const MDCellType& cell) const {
        using Vector3D = Vector<PosScalarType, Dim>;
        const size_t numMolecule = getNumMolecule();

        const MDCellType cellWithoutH(cell.getLattice(), cell.getPos().bottomRows(2 * numMolecule), cell.getMassVec());
        const ScalarType factor = 4 * epsilon;
        const ScalarType interMoleculeEnergy = LJModel.potentialEnergy(cellWithoutH) * factor;
        ScalarType intraMoleculeEnergy = 0;
        /* Intra molecule */ {
            Vector3D vecOH1, vecOH2;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                vecOH1 = cell.minDistVector(offset + i, 2 * i);
                vecOH2 = cell.minDistVector(offset + i, 2 * i + 1);
                const ScalarType r1 = vecOH1.norm();
                const ScalarType r2 = vecOH2.norm();
                const ScalarType elastic = square(arccos((vecOH1 * vecOH2) / (r1 * r2)) - ScalarType(equalTheta)) * (kTheta * 0.5);
                intraMoleculeEnergy += modifiedMorsePot(r1) + modifiedMorsePot(r2) + elastic;
            }
        }
        return interMoleculeEnergy + intraMoleculeEnergy;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::ewaldEnergy(const MDCellType& cell) const {
        const PositionMatrix chargePos = makeChargePos(cell);
        const size_t numMolecule = getNumMolecule();
        ScalarType selfCoulomb = 0;
        /* Spurious Coulomb Part */ {
            PeriodicCell<PosScalarType, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
            const size_t minIndexO = 2 * numMolecule;
            const size_t maxIndexO = minIndexO + numMolecule;
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                const size_t indexH1 = 2 * (i - minIndexO);
                const size_t indexH2 = indexH1 + 1;
                selfCoulomb += ScalarType(-charge * charge * 0.5) / chargeCell.minDistVector(i, indexH1).norm();
                selfCoulomb += ScalarType(-charge * charge * 0.5) / chargeCell.minDistVector(i, indexH2).norm();
                selfCoulomb += ScalarType(charge * charge * 0.25) / chargeCell.minDistVector(indexH1, indexH2).norm();
            }
        }
        return ewald.potentialEnergy(chargePos) - selfCoulomb;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::modifiedMorsePot(ScalarType r) {
        const ScalarType delta = (r - ScalarType(equalR)) * alphaR;
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        const ScalarType delta4 = square(delta2);
        return (delta2 - delta3 + delta4 * (7.0 / 12)) * Dr;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::modifiedMorseForce(ScalarType r) {
        const ScalarType delta = (r - ScalarType(equalR)) * alphaR;
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        return -(delta * 2 - delta2 * 3 + delta3 * (7.0 / 3)) * (Dr * alphaR);
    }
}
