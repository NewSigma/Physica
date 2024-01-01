/*
 * Copyright 2022-2024 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/ReshapedVector.h"
#include "Physica/Core/Physics/Molecular/WaterPolarTensor.h"
#include "Ewald/RSpaceEwald.h"
#include "LJModel.h"
#include "ForceModelImpl/AABModel.h"

namespace Physica::Core {
    template<class ScalarType, class EwaldType> class Q_TIP4P;

    namespace Internal {
        template<class ScalarType, class EwaldType>
        class Traits<Q_TIP4P<ScalarType, EwaldType>> {
        public:
            constexpr static bool IsPeriodBoundary = true;
            constexpr static bool IsLatticeDependent = true;
        };
    }
    /**
     * Reference:
     * [1] J. Chem. Phys. 131, 024501 (2009); https://doi.org/10.1063/1.3167790
     * [2] Jos Thijssen. Computational Physics[M]. London: Cambridge university press, 2013:205
     */
    template<class ScalarType, class EwaldType>
    class Q_TIP4P : public AABModel<ScalarType> {
        using This = Q_TIP4P<ScalarType, EwaldType>;
        using Base = AABModel<ScalarType>;
        using Base::Dim;
        using PlainScalar = typename ScalarType::PlainScalar;
        using REwaldType = typename Internal::Traits<EwaldType>::REwaldType;
        constexpr static bool IsSmallCell = Internal::Traits<REwaldType>::IsSmallCell;
        using LJModelType = LJModel<ScalarType, IsSmallCell>;
        using Vector3D = Vector<ScalarType, Dim>;
    public:
        using MDCellType = MDCell<ScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;

        constexpr static double epsilon = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(0.1852 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double epsilon4 = 4 * epsilon;
        constexpr static double lj_sigma = PhyConst<AU>::angstormToBohr(3.1589);
        constexpr static double charge = 1.1128;
        constexpr static double gamma = 0.73612;
        constexpr static double Dr = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(116.09 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double alphaR = PhyConst<AU>::bohrToAngstorm(2.287);
        constexpr static double equalR = PhyConst<AU>::angstormToBohr(0.9419);
        constexpr static double kTheta = PhyConst<AU>::eVToHartree(PhyConst<SI>::calorieToJoule(87.85 * 1000) / PhyConst<SI>::unitCharge) / PhyConst<SI>::avogadroNa;
        constexpr static double equalTheta = PhyConst<SI>::degreeToRadian(107.4);
    private:
        EwaldType ewald;
        LJModelType lj_model;
    public:
        Q_TIP4P(const MDCellType& refer_cell, PlainScalar cutoff, EwaldType ewald_);
        Q_TIP4P(const Q_TIP4P&) = default;
        Q_TIP4P(Q_TIP4P&&) noexcept = default;
        ~Q_TIP4P() = default;
        /* Operators */
        Q_TIP4P& operator=(Q_TIP4P model) noexcept;
        /* Operations */
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        [[nodiscard]] ScalarType potentialEnergy_unsort(const MDCellType& cell) const;

        template<class Executor> [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell);
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_unsort(const MDCellType& cell);
        template<class VectorType, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result);
        template<class Executor> [[nodiscard]] inline Vector<ScalarType> force_short(const MDCellType& cell);
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_short_unsort(const MDCellType& cell);
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell);
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_long_unsort(const MDCellType& cell);

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell);
        [[nodiscard]] LatticeMatrix virial_morse(const MDCellType& cell) const;
        template<class Executor, bool UseDynamicPolar>
        [[nodiscard]] PositionMatrix makeInducedDipole(const MDCellType& cell);
        void swap(Q_TIP4P& __restrict model) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return getNumParticle() / 3; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return ewald.getNumParticle(); }
        [[nodiscard]] const EwaldType& getEwald() const noexcept { return ewald; }
        [[nodiscard]] const typename MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Setters */
        inline void setLattice(LatticeMatrix lattice);
        /* Static members */
        [[nodiscard]] static PositionMatrix makePermanentDipole(const PeriodicCell<ScalarType, 3>& cell);
        inline static PermutationMatrix<ScalarType> sortPosition(MDCellType& cell);
    private:
        PositionMatrix makeChargePos(const MDCellType& cell) const;

        ScalarType potentialEnergyWithoutEwald(const MDCellType& cell) const;
        ScalarType ewaldEnergy(const MDCellType& cell) const;

        template<class Executor> void force_short_interMolecule(const MDCellType& cell, Vector<ScalarType>& shortForce) const;
        void force_short_intraMolecule(const MDCellType& cell, Vector<ScalarType>& shortForce) const;
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_short_PartialChargeRepr(const PositionMatrix& chargePos);
        template<class Executor> [[nodiscard]] inline Vector<ScalarType> force_long_PartialChargeRepr(const PositionMatrix& chargePos);
        template<class VectorType> void changeRepr(ContinuousVector<VectorType>& ewaldForce) const;

        [[nodiscard]] static ScalarType modifiedMorsePot(ScalarType r);
        [[nodiscard]] static ScalarType modifiedMorseForce(ScalarType r);
        [[nodiscard]] inline static MDCellType makeCellWithoutH(const MDCellType& original);
        [[nodiscard]] inline static bool isCellOrdered(const MDCellType& cell);
    };

    template<class ScalarType, class EwaldType>
    Q_TIP4P<ScalarType, EwaldType>::Q_TIP4P(const MDCellType& refer_cell, PlainScalar cutoff, EwaldType ewald_)
            : ewald(std::move(ewald_))
            , lj_model(lj_sigma, cutoff.getValue()) {
        assert(refer_cell.getNumParticle() % 3 == 0 && "[Error]: This is not a cell of water");
        const size_t numMolecule = refer_cell.getNumParticle() / 3;
        ewald = REwaldType(refer_cell.getLattice(), Base::makeCharges(numMolecule, charge * 0.5, -charge));
    }

    template<class ScalarType, class EwaldType>
    Q_TIP4P<ScalarType, EwaldType>& Q_TIP4P<ScalarType, EwaldType>::operator=(Q_TIP4P model) noexcept {
        swap(model);
        return *this;
    }

    template<class ScalarType, class EwaldType>
    ScalarType Q_TIP4P<ScalarType, EwaldType>::potentialEnergy(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Inconsistent atom number");
        assert(isCellOrdered(cell));
        return potentialEnergyWithoutEwald(cell) + ewaldEnergy(cell);
    }

    template<class ScalarType, class EwaldType>
    ScalarType Q_TIP4P<ScalarType, EwaldType>::potentialEnergy_unsort(const MDCellType& cell) const {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        return potentialEnergy(cell);
    }

    template<class ScalarType, class EwaldType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force(const MDCellType& cell) {
        Vector<ScalarType> result;
        forceAsync<Vector<ScalarType>, Executor>(cell, result);
        return result;
    }

    template<class ScalarType, class EwaldType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_unsort(const MDCellType& cell) {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const Vector<ScalarType> sort_f = force<Executor>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }

    template<class ScalarType, class EwaldType>
    template<class VectorType, class Executor>
    void Q_TIP4P<ScalarType, EwaldType>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) {
        assert(cell.getNumParticle() % 3 == 0);
        Vector<ScalarType> temp(result.getLength(), 0);
        auto future = Executor::schedule([this, &cell, &temp]() {
            force_short_interMolecule<Executor>(cell, temp);
            force_short_intraMolecule(cell, temp);
        });

        const auto chargePos = makeChargePos(cell);
        result = force_long_PartialChargeRepr<Executor>(chargePos);
        result += force_short_PartialChargeRepr<Executor>(chargePos);
        changeRepr(result);
        Executor::auto_wait(future);
        result += temp;
    }

    template<class ScalarType, class EwaldType>
    template<class Executor>
    inline Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_short(const MDCellType& cell) {
        assert(cell.getNumParticle() % 3 == 0);
        assert(isCellOrdered(cell));
        const size_t numMolecule = getNumMolecule();
        Vector<ScalarType> result(3 * numMolecule * Dim, 0);
        force_short_interMolecule<Executor>(cell, result);
        force_short_intraMolecule(cell, result);
        return result;
    }

    template<class ScalarType, class EwaldType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_short_unsort(const MDCellType& cell) {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const Vector<ScalarType> sort_f = force_short<Executor>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }

    template<class ScalarType, class EwaldType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_long(const MDCellType& cell) {
        auto result = force_long_PartialChargeRepr<Executor>(makeChargePos(cell));
        result += force_short_PartialChargeRepr<Executor>(makeChargePos(cell));
        changeRepr(result);
        return result;
    }

    template<class ScalarType, class EwaldType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_long_unsort(const MDCellType& cell) {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const Vector<ScalarType> sort_f = force_long<Executor>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }

    template<class ScalarType, class EwaldType>
    typename Q_TIP4P<ScalarType, EwaldType>::LatticeMatrix
    Q_TIP4P<ScalarType, EwaldType>::virial(const MDCellType& cell) {
        LatticeMatrix result(Dim, Dim, 0);
        const size_t numMolecule = getNumMolecule();
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = offset + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            Vector3D vecOH1 = cell.minDistVector(indexO, indexH1);
            Vector3D vecOH2 = cell.minDistVector(indexO, indexH2);
            const ScalarType r1 = vecOH1.norm();
            const ScalarType r2 = vecOH2.norm();
            /* Morse contribution */
            result += (vecOH1 * (modifiedMorseForce(r1) / r1)) * vecOH1.transpose();
            result += (vecOH2 * (modifiedMorseForce(r2) / r2)) * vecOH2.transpose();
            /* Angle contribution */
            LatticeMatrix temp = vecOH1 * vecOH2.transpose();
            LatticeMatrix angle = temp + temp.transpose();

            const ScalarType dot = vecOH1 * vecOH2;
            vecOH1 *= reciprocal(r1);
            vecOH2 *= reciprocal(r2);
            temp = vecOH1 * vecOH1.transpose() + vecOH2 * vecOH2.transpose();
            angle -= temp * dot;

            const ScalarType theta = arccos(dot / (r1 * r2));
            const ScalarType factor = (theta - ScalarType(equalTheta)) * (ScalarType(kTheta) / sqrt(square(r1 * r2) - square(dot)));
            result += angle * factor;
        }
        const PositionMatrix chargePos = makeChargePos(cell);
        /* Spurious Coulomb Part */ {
            PeriodicCell<ScalarType, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
            for (size_t i = 0; i < numMolecule; ++i) {
                const size_t indexO = offset + i;
                const size_t indexH1 = 2 * i;
                const size_t indexH2 = 2 * i + 1;
                const Vector3D vecOH1 = chargeCell.minDistVector(indexO, indexH1);
                const Vector3D vecOH2 = chargeCell.minDistVector(indexO, indexH2);
                const Vector3D vecH1H2 = vecOH2 - vecOH1;
                const ScalarType r1 = vecOH1.norm();
                const ScalarType r2 = vecOH2.norm();
                const ScalarType r12 = vecH1H2.norm();
                result -= (vecOH1 * (ScalarType(-charge * charge * 0.5) / (r1 * square(r1)))) * vecOH1.transpose();
                result -= (vecOH2 * (ScalarType(-charge * charge * 0.5) / (r2 * square(r2)))) * vecOH2.transpose();
                result -= (vecH1H2 * (ScalarType(charge * charge * 0.25) / (r12 * square(r12)))) * vecH1H2.transpose();
            }
        }
        result *= reciprocal(ewald.getVolume());

        result += ewald.virial(chargePos);
        result += lj_model.virial(makeCellWithoutH(cell)) * ScalarType(epsilon4);
        return result;
    }

    template<class ScalarType, class EwaldType>
    typename Q_TIP4P<ScalarType, EwaldType>::LatticeMatrix Q_TIP4P<ScalarType, EwaldType>::virial_morse(const MDCellType& cell) const {
        LatticeMatrix result(Dim, Dim, 0);
        const size_t numMolecule = getNumMolecule();
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = offset + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            Vector3D vecOH1 = cell.minDistVector(indexO, indexH1);
            Vector3D vecOH2 = cell.minDistVector(indexO, indexH2);
            const ScalarType r1 = vecOH1.norm();
            const ScalarType r2 = vecOH2.norm();
            /* Morse contribution */
            result += (vecOH1 * (modifiedMorseForce(r1) / r1)) * vecOH1.transpose();
            result += (vecOH2 * (modifiedMorseForce(r2) / r2)) * vecOH2.transpose();
        }
        result *= reciprocal(ewald.getVolume());
        return result;
    }

    template<class ScalarType, class EwaldType>
    template<class Executor, bool UseDynamicPolar>
    typename Q_TIP4P<ScalarType, EwaldType>::PositionMatrix
    Q_TIP4P<ScalarType, EwaldType>::makeInducedDipole(const MDCellType& cell) {
        assert(cell.getNumParticle() % 3 == 0 && "[Error]: This is not cell for water");
        const size_t numMolecule = cell.getNumParticle() / 3;
        const PositionMatrix chargePos = makeChargePos(cell);
        Vector<ScalarType> coulomb = force_short_PartialChargeRepr<Executor>(chargePos) + force_long_PartialChargeRepr<Executor>(chargePos);
        changeRepr(coulomb);
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

    template<class ScalarType, class EwaldType>
    void Q_TIP4P<ScalarType, EwaldType>::swap(Q_TIP4P& __restrict model) noexcept {
        assert(this != &model && "[Error]: Self swap is likely a bug");
        ewald.swap(model.ewald);
        lj_model.swap(model.lj_model);
    }

    template<class ScalarType, class EwaldType>
    inline void Q_TIP4P<ScalarType, EwaldType>::setLattice(LatticeMatrix lattice) {
        ewald.setLattice(std::move(lattice));
    }

    template<class ScalarType, class EwaldType>
    typename Q_TIP4P<ScalarType, EwaldType>::PositionMatrix
    Q_TIP4P<ScalarType, EwaldType>::makePermanentDipole(const PeriodicCell<ScalarType, 3>& cell) {
        assert(cell.getNumParticle() % 3 == 0 && "[Error]: This is not cell for water");
        const size_t numMolecule = cell.getNumParticle() / 3;
        PositionMatrix dipoles(numMolecule, 3);
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = 2 * numMolecule + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            auto dipole = dipoles.row(i);
            dipole = (cell.minDistVector(indexO, indexH1) + cell.minDistVector(indexO, indexH2)) * ScalarType(gamma * charge * 0.5);
        }
        return dipoles;
    }

    template<class ScalarType, class EwaldType>
    inline PermutationMatrix<ScalarType> Q_TIP4P<ScalarType, EwaldType>::sortPosition(MDCellType& cell) {
        return Base::sortPosition(cell, 1, 8);
    }

    template<class ScalarType, class EwaldType>
    typename Q_TIP4P<ScalarType, EwaldType>::PositionMatrix
    Q_TIP4P<ScalarType, EwaldType>::makeChargePos(const MDCellType& cell) const {
        assert(cell.getNumParticle() % 3 == 0);
        assert(isCellOrdered(cell));
        assert(ewald.getLattice() == cell.getLattice() && "[Error]: Lattice is not updated");

        PositionMatrix chargePos(cell.getPos().getRow(), 3);
        const auto& pos = cell.getPos();
        const size_t numMolecule = getNumMolecule();
        const size_t maxIndexH = 2 * numMolecule;
        auto chargePosH = chargePos.topRows(maxIndexH);
        chargePosH = pos.topRows(maxIndexH);

        const size_t minIndexO = maxIndexH;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType, 3> vecOH1, vecOH2;
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

    template<class ScalarType, class EwaldType>
    ScalarType Q_TIP4P<ScalarType, EwaldType>::potentialEnergyWithoutEwald(const MDCellType& cell) const {
        using Vector3D = Vector<ScalarType, Dim>;
        const size_t numMolecule = getNumMolecule();

        const ScalarType interMoleculeEnergy = lj_model.potentialEnergy(makeCellWithoutH(cell)) * ScalarType(epsilon4);
        ScalarType intraMoleculeEnergy = 0;
        /* Intra molecule */ {
            Vector3D vecOH1, vecOH2;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                vecOH1 = cell.minDistVector(offset + i, 2 * i);
                vecOH2 = cell.minDistVector(offset + i, 2 * i + 1);
                const ScalarType r1 = vecOH1.norm();
                const ScalarType r2 = vecOH2.norm();
                const ScalarType u = (vecOH1 * vecOH2) / (r1 * r2);
                const ScalarType elastic = square(arccos(u) - ScalarType(equalTheta)) * ScalarType(kTheta * 0.5);
                intraMoleculeEnergy += modifiedMorsePot(r1) + modifiedMorsePot(r2) + elastic;
            }
        }
        return interMoleculeEnergy + intraMoleculeEnergy;
    }

    template<class ScalarType, class EwaldType>
    ScalarType Q_TIP4P<ScalarType, EwaldType>::ewaldEnergy(const MDCellType& cell) const {
        assert(ewald.getLattice() == cell.getLattice() && "[Error]: Lattice is not updated");
        const PositionMatrix chargePos = makeChargePos(cell);
        const size_t numMolecule = getNumMolecule();
        ScalarType selfCoulomb = 0;
        /* Spurious Coulomb Part */ {
            PeriodicCell<ScalarType, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
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

    template<class ScalarType, class EwaldType>
    template<class Executor>
    void Q_TIP4P<ScalarType, EwaldType>::force_short_interMolecule(const MDCellType& cell, Vector<ScalarType>& shortForce) const {
        const size_t numMolecule = getNumMolecule();
        auto force_oxygen = shortForce.tail(2 * numMolecule * Dim);
        force_oxygen += lj_model.template force<Executor>(makeCellWithoutH(cell)) * ScalarType(epsilon4);
    }

    template<class ScalarType, class EwaldType>
    void Q_TIP4P<ScalarType, EwaldType>::force_short_intraMolecule(const MDCellType& cell, Vector<ScalarType>& shortForce) const {
        Vector3D vecOH1, vecOH2;
        Vector<ScalarType, 3> f;
        const size_t numMolecule = getNumMolecule();
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = offset + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            vecOH1 = cell.minDistVector(indexO, indexH1);
            vecOH2 = cell.minDistVector(indexO, indexH2);
            const ScalarType r1 = vecOH1.norm();
            const ScalarType r2 = vecOH2.norm();
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
            const ScalarType rep_r1 = reciprocal(vecOH1.norm());
            const ScalarType rep_r2 = reciprocal(vecOH2.norm());
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

    template<class ScalarType, class EwaldType>
    template<class Executor>
    Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_short_PartialChargeRepr(const PositionMatrix& chargePos) {
        Vector<ScalarType> coulomb = ewald.template force_short<Executor>(chargePos);
        PeriodicCell<ScalarType, 3> chargeCell(ewald.getLattice(), chargePos, MDCellType::Type::Cartesian);
        const size_t numMolecule = getNumMolecule();
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
            ScalarType r2 = vecMH1.squaredNorm();
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

    template<class ScalarType, class EwaldType>
    template<class Executor>
    inline Vector<ScalarType> Q_TIP4P<ScalarType, EwaldType>::force_long_PartialChargeRepr(const PositionMatrix& chargePos) {
        return ewald.template force_long<Executor>(chargePos);
    }
    /**
     * Change representation: from partial charge representation to HOH representation
     */
    template<class ScalarType, class EwaldType>
    template<class VectorType>
    void Q_TIP4P<ScalarType, EwaldType>::changeRepr(ContinuousVector<VectorType>& ewaldForce) const {
        const size_t numMolecule = getNumMolecule();
        const size_t minIndexO = 2 * numMolecule;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector<ScalarType, 3> f;
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            const size_t indexH1 = 2 * (i - minIndexO);
            const size_t indexH2 = indexH1 + 1;
            auto forceO = ewaldForce.segment(3 * i, 3 * i + 3);
            auto forceH1 = ewaldForce.segment(3 * indexH1, 3 * indexH1 + 3);
            auto forceH2 = ewaldForce.segment(3 * indexH2, 3 * indexH2 + 3);
            f = forceO * ScalarType((1 - gamma) * 0.5);
            forceH1 += f;
            forceH2 += f;
            forceO *= ScalarType(gamma);
        }
    }

    template<class ScalarType, class EwaldType>
    ScalarType Q_TIP4P<ScalarType, EwaldType>::modifiedMorsePot(ScalarType r) {
        const ScalarType delta = (r - ScalarType(equalR)) * ScalarType(alphaR);
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        const ScalarType delta4 = square(delta2);
        return (delta2 - delta3 + delta4 * ScalarType(7.0 / 12)) * ScalarType(Dr);
    }

    template<class ScalarType, class EwaldType>
    ScalarType Q_TIP4P<ScalarType, EwaldType>::modifiedMorseForce(ScalarType r) {
        const ScalarType delta = (r - ScalarType(equalR)) * ScalarType(alphaR);
        const ScalarType delta2 = square(delta);
        const ScalarType delta3 = delta2 * delta;
        return -(delta * ScalarType(2) - delta2 * ScalarType(3) + delta3 * ScalarType(7.0 / 3)) * ScalarType(Dr * alphaR);
    }

    template<class ScalarType, class EwaldType>
    inline typename Q_TIP4P<ScalarType, EwaldType>::MDCellType Q_TIP4P<ScalarType, EwaldType>::makeCellWithoutH(
            const MDCellType& original) {
        const size_t numMolecule = original.getNumParticle() / 3;
        return MDCellType(original.getLattice(), original.getPos().bottomRows(2 * numMolecule), original.getMassVec());
    }

    template<class ScalarType, class EwaldType>
    inline bool Q_TIP4P<ScalarType, EwaldType>::isCellOrdered(const MDCellType& cell) {
        return Base::isCellOrdered(cell, 1, 8);
    }
}
