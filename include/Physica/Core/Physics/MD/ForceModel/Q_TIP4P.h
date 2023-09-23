/*
 * Copyright 2022-2023 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermutationMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/ReshapedVector.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/Molecular/WaterPolarTensor.h"
#include "LJModel.h"
#include "Ewald.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType> class Q_TIP4P;
    /**
     * Reference:
     * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
     * [2] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:205
     */
    template<class ScalarType, class PosScalarType>
    class Q_TIP4P {
        constexpr static unsigned int Dim = 3;
        using This = Q_TIP4P<ScalarType, PosScalarType>;
        using EwaldType = Ewald<ScalarType, PosScalarType>;
        using LJModelType = LJModel<ScalarType, PosScalarType>;
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;

        constexpr static double epsilon = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(0.1852 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double charge = 1.1128;
        constexpr static double gamma = 0.73612;
        constexpr static double Dr = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(116.09 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double alphaR = PhyConst<AU>::bohrToAngstorm(2.287);
        constexpr static double equalR = PhyConst<AU>::angstormToBohr(0.9419);
        constexpr static double kTheta = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(87.85 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double equalTheta = PhyConst<SI>::degreeToRadian(107.4);

        constexpr static double lj_sigma = PhyConst<AU>::angstormToBohr(3.1589);
    private:
        size_t numMolecule;
        EwaldType ewald;
        LJModelType lj_model;
    public:
        Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_);
        Q_TIP4P(const Q_TIP4P&) = default;
        Q_TIP4P(Q_TIP4P&&) noexcept = default;
        ~Q_TIP4P() = default;
        /* Operators */
        Q_TIP4P& operator=(Q_TIP4P model) noexcept;
        /* Operations */
        template<class Executor, bool IsSmallCell = false> [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class Executor, bool IsSmallCell = false> [[nodiscard]] Vector<ScalarType> force_unsort(const MDCellType& cell) const;
        template<class VectorType, class Executor, bool IsSmallCell>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor, bool IsSmallCell = false> [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const;
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const;

        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        [[nodiscard]] ScalarType potentialEnergy_unsort(const MDCellType& cell) const;

        template<class Executor, bool UseDynamicPolar>
        [[nodiscard]] PositionMatrix makeInducedDipole(const MDCellType& cell) const;
        void swap(Q_TIP4P& model) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return numMolecule; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return getNumMolecule() * 3; }
        [[nodiscard]] const typename MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Setters */
        inline void updateLattice(LatticeMatrix lattice);
        /* Static members */
        [[nodiscard]] static PositionMatrix makePermanentDipole(const PeriodicCell<PosScalarType, 3>& cell);
        static PermutationMatrix<PosScalarType> sortPosition(MDCellType& cell);
    private:
        Vector<ScalarType> makeCharges() const;
        PositionMatrix makeChargePos(const MDCellType& cell) const;
        void force_short_intraMolecule(const MDCellType& cell, Vector<ScalarType>& shortForce) const;
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_long_PartialChargeRepr(const MDCellType& cell) const;

        ScalarType potentialEnergyWithoutEwald(const MDCellType& cell) const;
        ScalarType ewaldEnergy(const MDCellType& cell) const;

        static ScalarType modifiedMorsePot(ScalarType r);
        static ScalarType modifiedMorseForce(ScalarType r);
        static bool isCellOrdered(const MDCellType& cell);
    };

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>::Q_TIP4P(const MDCellType& refer_cell, ScalarType cutoff_)
            : numMolecule(refer_cell.getNumParticle() / 3)
            , lj_model(lj_sigma, std::move(cutoff_)) {
        assert(refer_cell.getNumParticle() % 3 == 0);
        ewald = EwaldType(refer_cell.getLattice(), makeCharges());
    }

    template<class ScalarType, class PosScalarType>
    Q_TIP4P<ScalarType, PosScalarType>& Q_TIP4P<ScalarType, PosScalarType>::operator=(Q_TIP4P model) noexcept {
        swap(model);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force(const MDCellType& cell) const {
        Vector<ScalarType> result;
        forceAsync<Vector<ScalarType>, Executor, IsSmallCell>(cell, result);
        return result;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force_unsort(const MDCellType& cell) const {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const Vector<ScalarType> sort_f = force<Executor, IsSmallCell>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }

    template<class ScalarType, class PosScalarType>
    template<class VectorType, class Executor, bool IsSmallCell>
    void Q_TIP4P<ScalarType, PosScalarType>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        static_assert(!Internal::Traits<Executor>::isCudaEnabled, "[Error]: Cuda is not supported");
        assert(cell.getNumParticle() % 3 == 0);
        auto future = Executor::schedule([this, &cell, &result]() {
            result = force_short<Executor, IsSmallCell>(cell);
        });

        const Vector<ScalarType> coulomb = force_long<Executor>(cell);
        Executor::auto_wait(future);
        result += coulomb;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force_short(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0);
        assert(isCellOrdered(cell));
        Vector<ScalarType> shortForce(3 * numMolecule * Dim, 0);
        /* LJ */ {
            const MDCellType cellWithoutH(cell.getLattice(), cell.getPos().bottomRows(2 * numMolecule), cell.getMassVec());
            const ScalarType factor = ScalarType(24 * epsilon / lj_sigma);
            auto force = shortForce.tail(2 * numMolecule * Dim);
            force = lj_model.template force<Executor, IsSmallCell>(cellWithoutH) * factor;
        }
        force_short_intraMolecule(cell, shortForce);
        return shortForce;
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force_long(const MDCellType& cell) const {
        Vector<ScalarType> coulomb = force_long_PartialChargeRepr<Executor>(cell);
        const size_t minIndexO = 2 * numMolecule;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType, 3> f;
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
        return coulomb;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergy(const MDCellType& cell) const {
        assert(isCellOrdered(cell));
        return potentialEnergyWithoutEwald(cell) + ewaldEnergy(cell);
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergy_unsort(const MDCellType& cell) const {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        return potentialEnergy(cell);
    }

    template<class ScalarType, class PosScalarType>
    template<class Executor, bool UseDynamicPolar>
    typename Q_TIP4P<ScalarType, PosScalarType>::PositionMatrix
    Q_TIP4P<ScalarType, PosScalarType>::makeInducedDipole(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0 && "[Error]: This is not cell for water");
        const size_t numMolecule = cell.getNumParticle() / 3;
        const Vector<ScalarType> coulomb = force_long_PartialChargeRepr<Executor>(cell);
        const ScalarType repCharge = reciprocal(ScalarType(charge));

        PositionMatrix dipoles(numMolecule, 3);
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = 2 * numMolecule + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            const auto posOH1 = cell.minDistVector(indexO, indexH1);
            const auto posOH2 = cell.minDistVector(indexO, indexH2);
            const WaterPolarTensor<ScalarType, UseDynamicPolar> polar(posOH1, posOH2);
            const Vector<ScalarType, 3> electricField = coulomb.segment(indexO * Dim, (indexO + 1) * Dim) * repCharge;
            auto dipole = dipoles.row(i);
            dipole = polar * electricField;
        }
        return dipoles;
    }

    template<class ScalarType, class PosScalarType>
    void Q_TIP4P<ScalarType, PosScalarType>::swap(Q_TIP4P& model) noexcept {
        assert(this != &model && "[Error]: Self swap is likely a bug");
        std::swap(numMolecule, model.numMolecule);
        ewald.swap(model.ewald);
        lj_model.swap(model.lj_model);
    }

    template<class ScalarType, class PosScalarType>
    inline void Q_TIP4P<ScalarType, PosScalarType>::updateLattice(LatticeMatrix lattice) {
        ewald.setLattice(std::move(lattice));
    }

    template<class ScalarType, class PosScalarType>
    typename Q_TIP4P<ScalarType, PosScalarType>::PositionMatrix
    Q_TIP4P<ScalarType, PosScalarType>::makePermanentDipole(const PeriodicCell<PosScalarType, 3>& cell) {
        assert(cell.getNumParticle() % 3 == 0 && "[Error]: This is not cell for water");
        const size_t numMolecule = cell.getNumParticle() / 3;
        PositionMatrix dipoles(numMolecule, 3);
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = 2 * numMolecule + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            auto dipole = dipoles.row(i);
            dipole = (cell.minDistVector(indexO, indexH1) + cell.minDistVector(indexO, indexH2)) * ScalarType(-gamma * charge);
        }
        return dipoles;
    }

    template<class ScalarType, class PosScalarType>
    PermutationMatrix<PosScalarType> Q_TIP4P<ScalarType, PosScalarType>::sortPosition(MDCellType& cell) {
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
                const bool isHydrogen = cell.getMass(i) == PhyConst<AU>::atomMass(1);
                const size_t index = isHydrogen ? indexH : indexO;
                new_pos.row(index) = source.row(i);
                new_mass[i] = i < numH ? PhyConst<AU>::atomMass(1) : PhyConst<AU>::atomMass(8);
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
                orderStage2[2 * i] = orderStage1[indexH1];

                auto posH2 = new_pos.row(2 * i + 1);
                posH2 = source.row(indexH2).asVector();
                orderStage2[2 * i + 1] = orderStage1[indexH2];
                orderStage2[indexO] = orderStage1[indexO];
            }
            cell = MDCellType(cell.getLattice(), std::move(new_pos), std::move(new_mass));
        }
        return PermutationMatrix<PosScalarType>(std::move(orderStage2));
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
        
        for (size_t i = minIndexO; i < maxIndexO; ++i)
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
    void Q_TIP4P<ScalarType, PosScalarType>::force_short_intraMolecule(const MDCellType& cell, Vector<ScalarType>& shortForce) const {
        using Vector3D = Vector<PosScalarType, Dim>;
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

    template<class ScalarType, class PosScalarType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, PosScalarType>::force_long_PartialChargeRepr(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0);
        assert(isCellOrdered(cell));
        assert(ewald.getLattice() == cell.getLattice() && "[Error]: Lattice is not updated");
        using Vector3D = Vector<PosScalarType, Dim>;

        const PositionMatrix chargePos = makeChargePos(cell);
        Vector<ScalarType> coulomb = ewald.template force<Executor>(chargePos);
        PeriodicCell<PosScalarType, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
        const size_t minIndexO = 2 * numMolecule;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType, 3> f;
        Vector3D vecMH1, vecMH2, vecH1H2;
        /* Remove self interaction */
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
        return coulomb;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType Q_TIP4P<ScalarType, PosScalarType>::potentialEnergyWithoutEwald(const MDCellType& cell) const {
        using Vector3D = Vector<PosScalarType, Dim>;
        const size_t numMolecule = getNumMolecule();

        const MDCellType cellWithoutH(cell.getLattice(), cell.getPos().bottomRows(2 * numMolecule), cell.getMassVec());
        const ScalarType factor = 4 * epsilon;
        const ScalarType interMoleculeEnergy = lj_model.potentialEnergy(cellWithoutH) * factor;
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
        assert(ewald.getLattice() == cell.getLattice() && "[Error]: Lattice is not updated");
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

    template<class ScalarType, class PosScalarType>
    bool Q_TIP4P<ScalarType, PosScalarType>::isCellOrdered(const MDCellType& cell) {
        const size_t numParticle= cell.getNumParticle();
        const size_t maxIndexH = 2 * numParticle / 3;
        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numParticle / 3;
        for (size_t i = 0; i < maxIndexH; ++i) {
            const bool isHydrogen = cell.getMass(i) == PhyConst<AU>::atomMass(1);
            if (!isHydrogen)
                return false;
        }
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            const bool isOxygen = cell.getMass(i) == PhyConst<AU>::atomMass(8);
            if (!isOxygen)
                return false;
        }
        return true;
    }
}
