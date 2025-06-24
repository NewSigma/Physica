/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Physics/Molecular/WaterPolarTensor.h"
#include "ForceModelImpl/AABModel.h"
#include "LJModel.h"

namespace Physica {
    /**
     * Reference:
     * [1] J. Chem. Phys. 131, 024501 (2009); https://doi.org/10.1063/1.3167790
     * [2] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:205
     */
    template<Scalar T, class EwaldType>
    class Q_TIP4P : public AABModel<T> {
        using This = Q_TIP4P<T, EwaldType>;
        using Base = AABModel<T>;
        using Base::Dim;
        using Tv = T::ValueType;
        using REwaldType = Traits<EwaldType>::REwaldType;
        using BornChargeArray = REwaldType::BornChargeArray;
        constexpr static bool IsSmallCell = Traits<REwaldType>::IsSmallCell;
        using LJModelType = LJModel<T, IsSmallCell>;

        static_assert(Dim == 3, "[Error]: Not implemented");
    public:
        using MDCellType = MDCell<T, Dim>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;

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
        Q_TIP4P(const MDCellType& refer_cell, Tv cutoff, EwaldType ewald_);
        Q_TIP4P(const This&) = default;
        Q_TIP4P(This&&) noexcept = default;
        ~Q_TIP4P() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T potentialV(const MDCellType& cell) const;
        [[nodiscard]] T potentialV_unsort(const MDCellType& cell) const;

        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force(const MDCellType& cell);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_unsort(const MDCellType& cell);
        template<Vector V, ExecutePolicy P>
        void forceAsync(const MDCellType& cell, ContinuousVector<V>& result);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_short(const MDCellType& cell);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_short_unsort(const MDCellType& cell);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_long(const MDCellType& cell);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_long_unsort(const MDCellType& cell);

        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_contract(const MDCellType& cell);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_uncontract(const MDCellType& cell);

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell);
        [[nodiscard]] LatticeMatrix virial_morse(const MDCellType& cell) const;
        [[nodiscard]] BornChargeArray calcBornCharge() const;
        template<ExecutePolicy P, bool UseDynamicPolar>
        [[nodiscard]] PositionMatrix makeInducedDipole(const MDCellType& cell);
        void swap(This& __restrict model) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return getNumParticle() / 3; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return ewald.getNumParticle(); }
        [[nodiscard]] const EwaldType& getEwald() const noexcept { return ewald; }
        [[nodiscard]] const MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Setters */
        void setLattice(LatticeMatrix lattice);
        /* Static members */
        [[nodiscard]] static PositionMatrix makePermanentDipole(const PeriodicCell<T, 3>& cell);
        static PermMatrix<T> sortPosition(MDCellType& cell);
    private:
        PositionMatrix makeChargePos(const MDCellType& cell) const;

        T potentialVWithoutEwald(const MDCellType& cell) const;
        T ewaldEnergy(const MDCellType& cell) const;

        template<ExecutePolicy P> void force_short_interMolecule(const MDCellType& cell, VectorND<T>& shortForce) const;
        void force_short_intraMolecule(const MDCellType& cell, VectorND<T>& shortForce) const;
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_short_PartialChargeRepr(const PositionMatrix& chargePos);
        template<ExecutePolicy P> [[nodiscard]] VectorND<T> force_long_PartialChargeRepr(const PositionMatrix& chargePos);
        template<Vector V> void changeRepr(ContinuousVector<V>& ewaldForce) const;

        [[nodiscard]] static T modifiedMorsePot(T r);
        [[nodiscard]] static T modifiedMorseForce(T r);
        [[nodiscard]] static MDCellType makeCellWithoutH(const MDCellType& original);
        [[nodiscard]] static bool isCellOrdered(const MDCellType& cell);
    };

    template<Scalar T, class EwaldType>
    Q_TIP4P<T, EwaldType>::Q_TIP4P(const MDCellType& refer_cell, Tv cutoff, EwaldType ewald_)
            : ewald(std::move(ewald_))
            , lj_model(lj_sigma, cutoff.value()) {
        assert(refer_cell.getNumParticle() % 3 == 0 && "[Error]: This is not a cell of water");
        const size_t numMolecule = refer_cell.getNumParticle() / 3;
        ewald = REwaldType(refer_cell.getLattice(), Base::makeCharges(numMolecule, charge * 0.5, -charge));
    }

    template<Scalar T, class EwaldType>
    T Q_TIP4P<T, EwaldType>::potentialV(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Inconsistent atom number");
        assert(isCellOrdered(cell));
        return potentialVWithoutEwald(cell) + ewaldEnergy(cell);
    }

    template<Scalar T, class EwaldType>
    T Q_TIP4P<T, EwaldType>::potentialV_unsort(const MDCellType& cell) const {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        return potentialV(cell);
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force(const MDCellType& cell) {
        VectorND<T> result;
        forceAsync<VectorND<T>, P>(cell, result);
        return result;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_unsort(const MDCellType& cell) {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const VectorND<T> sort_f = force<P>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }

    template<Scalar T, class EwaldType>
    template<Vector V, ExecutePolicy P>
    void Q_TIP4P<T, EwaldType>::forceAsync(const MDCellType& cell, ContinuousVector<V>& result) {
        assert(cell.getNumParticle() % 3 == 0);
        VectorND<T> temp(getNumParticle() * 3, 0);
        auto task = schedule<P>([this, &cell, &temp]() {
            force_short_interMolecule<P>(cell, temp);
            force_short_intraMolecule(cell, temp);
        });

        const auto chargePos = makeChargePos(cell);
        result = force_long_PartialChargeRepr<P>(chargePos);
        result += force_short_PartialChargeRepr<P>(chargePos);
        changeRepr(result);
        task.wait();
        result += temp;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_short(const MDCellType& cell) {
        assert(cell.getNumParticle() % 3 == 0);
        assert(isCellOrdered(cell));
        VectorND<T> result = force_short_PartialChargeRepr<P>(makeChargePos(cell));
        changeRepr(result);
        force_short_interMolecule<P>(cell, result);
        force_short_intraMolecule(cell, result);
        return result;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_short_unsort(const MDCellType& cell) {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const VectorND<T> sort_f = force_short<P>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_long(const MDCellType& cell) {
        auto result = force_long_PartialChargeRepr<P>(makeChargePos(cell));
        changeRepr(result);
        return result;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_long_unsort(const MDCellType& cell) {
        MDCellType copy = cell;
        const auto permute = sortPosition(copy);
        const VectorND<T> sort_f = force_long<P>(copy);
        const PositionMatrix unsort_f = permute.inverse() * sort_f.reshape(cell.getPos());
        return unsort_f.flatten();
    }
    /**
     * Contract strategy as [1] used
     * 
     * Reference:
     * [1] J. Chem. Phys. 129, 024105 (2008); https://doi.org/10.1063/1.2953308
     */
    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_contract(const MDCellType& cell) {
        auto result = force_long_PartialChargeRepr<P>(makeChargePos(cell));
        result += force_short_PartialChargeRepr<P>(makeChargePos(cell));
        changeRepr(result);
        force_short_interMolecule<P>(cell, result);
        return result;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_uncontract(const MDCellType& cell) {
        VectorND<T> result(3 * getNumMolecule() * Dim, 0);
        force_short_intraMolecule(cell, result);
        return result;
    }

    template<Scalar T, class EwaldType>
    auto Q_TIP4P<T, EwaldType>::virial(const MDCellType& cell) -> LatticeMatrix {
        LatticeMatrix result(Dim, Dim, 0);
        const size_t numMolecule = getNumMolecule();
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = offset + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            Vector3D<T> vecOH1 = cell.minDistVector(indexO, indexH1);
            Vector3D<T> vecOH2 = cell.minDistVector(indexO, indexH2);
            const T r1 = vecOH1.norm();
            const T r2 = vecOH2.norm();
            /* Morse contribution */
            result += (vecOH1 * (modifiedMorseForce(r1) / r1)) * vecOH1.transpose();
            result += (vecOH2 * (modifiedMorseForce(r2) / r2)) * vecOH2.transpose();
            /* Angle contribution */
            LatticeMatrix temp = vecOH1 * vecOH2.transpose();
            LatticeMatrix angle = temp + temp.transpose();

            const T dot = vecOH1 * vecOH2;
            vecOH1 *= reciprocal(r1);
            vecOH2 *= reciprocal(r2);
            temp = vecOH1 * vecOH1.transpose() + vecOH2 * vecOH2.transpose();
            angle -= temp * dot;

            const T theta = arccos(dot / (r1 * r2));
            const T factor = (theta - T(equalTheta)) * (T(kTheta) / sqrt(square(r1 * r2) - square(dot)));
            result += angle * factor;
        }
        const PositionMatrix chargePos = makeChargePos(cell);
        /* Spurious Coulomb Part */ {
            PeriodicCell<T, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
            for (size_t i = 0; i < numMolecule; ++i) {
                const size_t indexO = offset + i;
                const size_t indexH1 = 2 * i;
                const size_t indexH2 = 2 * i + 1;
                const Vector3D<T> vecOH1 = chargeCell.minDistVector(indexO, indexH1);
                const Vector3D<T> vecOH2 = chargeCell.minDistVector(indexO, indexH2);
                const Vector3D<T> vecH1H2 = vecOH2 - vecOH1;
                const T r1 = vecOH1.norm();
                const T r2 = vecOH2.norm();
                const T r12 = vecH1H2.norm();
                result -= (vecOH1 * (T(-charge * charge * 0.5) / (r1 * square(r1)))) * vecOH1.transpose();
                result -= (vecOH2 * (T(-charge * charge * 0.5) / (r2 * square(r2)))) * vecOH2.transpose();
                result -= (vecH1H2 * (T(charge * charge * 0.25) / (r12 * square(r12)))) * vecH1H2.transpose();
            }
        }
        result *= reciprocal(ewald.getVolume());

        result += ewald.virial(chargePos);
        result += lj_model.virial(makeCellWithoutH(cell)) * T(epsilon4);
        return result;
    }

    template<Scalar T, class EwaldType>
    auto Q_TIP4P<T, EwaldType>::virial_morse(const MDCellType& cell) const -> LatticeMatrix {
        LatticeMatrix result(Dim, Dim, 0);
        const size_t numMolecule = getNumMolecule();
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = offset + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            Vector3D<T> vecOH1 = cell.minDistVector(indexO, indexH1);
            Vector3D<T> vecOH2 = cell.minDistVector(indexO, indexH2);
            const T r1 = vecOH1.norm();
            const T r2 = vecOH2.norm();
            /* Morse contribution */
            result += (vecOH1 * (modifiedMorseForce(r1) / r1)) * vecOH1.transpose();
            result += (vecOH2 * (modifiedMorseForce(r2) / r2)) * vecOH2.transpose();
        }
        result *= reciprocal(ewald.getVolume());
        return result;
    }

    template<Scalar T, class EwaldType>
    auto Q_TIP4P<T, EwaldType>::calcBornCharge() const -> BornChargeArray {
        BornChargeArray result = ewald.calcBornCharge();
        for (size_t i = 0; i < result.getLength(); ++i) {
            auto diag = result[i].diag();
            diag *= T(gamma);
        }
        return result;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P, bool UseDynamicPolar>
    auto Q_TIP4P<T, EwaldType>::makeInducedDipole(const MDCellType& cell) -> PositionMatrix {
        assert(cell.getNumParticle() % 3 == 0 && "[Error]: This is not cell for water");
        const size_t numMolecule = cell.getNumParticle() / 3;
        const PositionMatrix chargePos = makeChargePos(cell);
        VectorND<T> coulomb = force_short_PartialChargeRepr<P>(chargePos) + force_long_PartialChargeRepr<P>(chargePos);
        changeRepr(coulomb);
        const T repCharge = reciprocal(T(charge));

        PositionMatrix dipoles(numMolecule, 3);
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = 2 * numMolecule + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            const auto posOH1 = cell.minDistVector(indexO, indexH1);
            const auto posOH2 = cell.minDistVector(indexO, indexH2);
            const WaterPolarTensor<T, UseDynamicPolar> polar(posOH1, posOH2);
            const Vector3D<T> electricField = coulomb.segment(indexO * Dim, (indexO + 1) * Dim) * repCharge;
            auto dipole = dipoles.row(i);
            dipole = polar * electricField;
        }
        return dipoles;
    }

    template<Scalar T, class EwaldType>
    void Q_TIP4P<T, EwaldType>::swap(This& __restrict model) noexcept {
        assert(this != &model && "[Error]: Self swap is likely a bug");
        ewald.swap(model.ewald);
        lj_model.swap(model.lj_model);
    }

    template<Scalar T, class EwaldType>
    void Q_TIP4P<T, EwaldType>::setLattice(LatticeMatrix lattice) {
        ewald.setLattice(std::move(lattice));
    }

    template<Scalar T, class EwaldType>
    auto Q_TIP4P<T, EwaldType>::makePermanentDipole(const PeriodicCell<T, 3>& cell) -> PositionMatrix {
        assert(cell.getNumParticle() % 3 == 0 && "[Error]: This is not cell for water");
        const size_t numMolecule = cell.getNumParticle() / 3;
        PositionMatrix dipoles(numMolecule, 3);
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = 2 * numMolecule + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            auto dipole = dipoles.row(i);
            dipole = (cell.minDistVector(indexO, indexH1) + cell.minDistVector(indexO, indexH2)) * T(gamma * charge * 0.5);
        }
        return dipoles;
    }

    template<Scalar T, class EwaldType>
    PermMatrix<T> Q_TIP4P<T, EwaldType>::sortPosition(MDCellType& cell) {
        return Base::sortPosition(cell, 1, 8);
    }

    template<Scalar T, class EwaldType>
    auto Q_TIP4P<T, EwaldType>::makeChargePos(const MDCellType& cell) const -> PositionMatrix {
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
        Vector3D<T> vecOH1, vecOH2;
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            auto chargePosO = chargePos.row(i);
            auto posO = pos.row(i);
            const size_t indexH1 = 2 * (i - minIndexO);
            const size_t indexH2 = indexH1 + 1;
            vecOH1 = cell.minDistVector(i, indexH1);
            vecOH2 = cell.minDistVector(i, indexH2);
            chargePosO = posO + (vecOH1 + vecOH2) * T((1 - gamma) * 0.5);
        }
        cell.normalizePos(chargePos);
        return chargePos;
    }

    template<Scalar T, class EwaldType>
    T Q_TIP4P<T, EwaldType>::potentialVWithoutEwald(const MDCellType& cell) const {
        const size_t numMolecule = getNumMolecule();

        const T interMoleculeEnergy = lj_model.potentialV(makeCellWithoutH(cell)) * T(epsilon4);
        T intraMoleculeEnergy = 0;
        /* Intra molecule */ {
            Vector3D<T> vecOH1, vecOH2;
            const size_t offset = 2 * numMolecule;
            for (size_t i = 0; i < numMolecule; ++i) {
                vecOH1 = cell.minDistVector(offset + i, 2 * i);
                vecOH2 = cell.minDistVector(offset + i, 2 * i + 1);
                const T r1 = vecOH1.norm();
                const T r2 = vecOH2.norm();
                const T u = (vecOH1 * vecOH2) / (r1 * r2);
                const T elastic = square(arccos(u) - T(equalTheta)) * T(kTheta * 0.5);
                intraMoleculeEnergy += modifiedMorsePot(r1) + modifiedMorsePot(r2) + elastic;
            }
        }
        return interMoleculeEnergy + intraMoleculeEnergy;
    }

    template<Scalar T, class EwaldType>
    T Q_TIP4P<T, EwaldType>::ewaldEnergy(const MDCellType& cell) const {
        assert(ewald.getLattice() == cell.getLattice() && "[Error]: Lattice is not updated");
        const PositionMatrix chargePos = makeChargePos(cell);
        const size_t numMolecule = getNumMolecule();
        T selfCoulomb = 0;
        /* Spurious Coulomb Part */ {
            PeriodicCell<T, 3> chargeCell(cell.getLattice(), chargePos, cell.getType());
            const size_t minIndexO = 2 * numMolecule;
            const size_t maxIndexO = minIndexO + numMolecule;
            for (size_t i = minIndexO; i < maxIndexO; ++i) {
                const size_t indexH1 = 2 * (i - minIndexO);
                const size_t indexH2 = indexH1 + 1;
                selfCoulomb += T(-charge * charge * 0.5) / chargeCell.minDistVector(i, indexH1).norm();
                selfCoulomb += T(-charge * charge * 0.5) / chargeCell.minDistVector(i, indexH2).norm();
                selfCoulomb += T(charge * charge * 0.25) / chargeCell.minDistVector(indexH1, indexH2).norm();
            }
        }
        return ewald.potentialV(chargePos) - selfCoulomb;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    void Q_TIP4P<T, EwaldType>::force_short_interMolecule(const MDCellType& cell, VectorND<T>& shortForce) const {
        const size_t numMolecule = getNumMolecule();
        auto force_oxygen = shortForce.tail(2 * numMolecule * Dim);
        force_oxygen += lj_model.template force<P>(makeCellWithoutH(cell)) * T(epsilon4);
    }

    template<Scalar T, class EwaldType>
    void Q_TIP4P<T, EwaldType>::force_short_intraMolecule(const MDCellType& cell, VectorND<T>& shortForce) const {
        Vector3D<T> vecOH1, vecOH2;
        Vector3D<T> f;
        const size_t numMolecule = getNumMolecule();
        const size_t offset = 2 * numMolecule;
        for (size_t i = 0; i < numMolecule; ++i) {
            const size_t indexO = offset + i;
            const size_t indexH1 = 2 * i;
            const size_t indexH2 = 2 * i + 1;
            vecOH1 = cell.minDistVector(indexO, indexH1);
            vecOH2 = cell.minDistVector(indexO, indexH2);
            const T r1 = vecOH1.norm();
            const T r2 = vecOH2.norm();
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

        Vector3D<T> force1, force2;
        for (size_t i = 0; i < numMolecule; ++i) {
            vecOH1 = cell.minDistVector(offset + i, 2 * i);
            vecOH2 = cell.minDistVector(offset + i, 2 * i + 1);
            const T rep_r1 = reciprocal(vecOH1.norm());
            const T rep_r2 = reciprocal(vecOH2.norm());
            const T u = (vecOH1 * vecOH2) * (rep_r1 * rep_r2);
            const T temp = reciprocal(sqrt(T(1) - square(u)));
            const T factor1 = square(rep_r1) * u * temp;
            const T factor2 = square(rep_r2) * u * temp;
            const T factor3 = rep_r1 * rep_r2 * temp;

            force1 = vecOH1 * factor1 - vecOH2 * factor3;
            force2 = vecOH2 * factor2 - vecOH1 * factor3;
            const T theta = arccos(u);
            const T factor = (T(equalTheta) - theta) * kTheta;
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

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_short_PartialChargeRepr(const PositionMatrix& chargePos) {
        VectorND<T> coulomb = ewald.template force_short<P>(chargePos);
        PeriodicCell<T, 3> chargeCell(ewald.getLattice(), chargePos, MDCellType::Type::Cartesian);
        const size_t numMolecule = getNumMolecule();
        const size_t minIndexO = 2 * numMolecule;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector3D<T> f;
        Vector3D<T> vecMH1, vecMH2, vecH1H2;
        /* Remove self interaction */
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            const size_t indexH1 = 2 * (i - minIndexO);
            const size_t indexH2 = indexH1 + 1;
            auto forceM = coulomb.segment(3 * i, 3 * i + 3);
            auto forceH1 = coulomb.segment(3 * indexH1, 3 * indexH1 + 3);
            auto forceH2 = coulomb.segment(3 * indexH2, 3 * indexH2 + 3);

            vecMH1 = chargeCell.minDistVector(i, indexH1);
            T r2 = vecMH1.squaredNorm();
            f = vecMH1 * (T(-charge * charge * 0.5) / (r2 * sqrt(r2)));
            forceM += f;
            forceH1 -= f;

            vecMH2 = chargeCell.minDistVector(i, indexH2);
            r2 = vecMH2.squaredNorm();
            f = vecMH2 * (T(-charge * charge * 0.5) / (r2 * sqrt(r2)));
            forceM += f;
            forceH2 -= f;

            vecH1H2 = chargeCell.minDistVector(indexH1, indexH2);
            r2 = vecH1H2.squaredNorm();
            f = vecH1H2 * (T(charge * charge * 0.25) / (r2 * sqrt(r2)));
            forceH1 += f;
            forceH2 -= f;
        }
        return coulomb;
    }

    template<Scalar T, class EwaldType>
    template<ExecutePolicy P>
    VectorND<T> Q_TIP4P<T, EwaldType>::force_long_PartialChargeRepr(const PositionMatrix& chargePos) {
        return ewald.template force_long<P>(chargePos);
    }
    /**
     * Change representation: from partial charge representation to HOH representation
     */
    template<Scalar T, class EwaldType>
    template<Vector V>
    void Q_TIP4P<T, EwaldType>::changeRepr(ContinuousVector<V>& ewaldForce) const {
        const size_t numMolecule = getNumMolecule();
        const size_t minIndexO = 2 * numMolecule;
        const size_t maxIndexO = minIndexO + numMolecule;
        Vector3D<T> f;
        for (size_t i = minIndexO; i < maxIndexO; ++i) {
            const size_t indexH1 = 2 * (i - minIndexO);
            const size_t indexH2 = indexH1 + 1;
            auto forceO = ewaldForce.segment(3 * i, 3 * i + 3);
            auto forceH1 = ewaldForce.segment(3 * indexH1, 3 * indexH1 + 3);
            auto forceH2 = ewaldForce.segment(3 * indexH2, 3 * indexH2 + 3);
            f = forceO * T((1 - gamma) * 0.5);
            forceH1 += f;
            forceH2 += f;
            forceO *= T(gamma);
        }
    }

    template<Scalar T, class EwaldType>
    T Q_TIP4P<T, EwaldType>::modifiedMorsePot(T r) {
        const T delta = (r - T(equalR)) * T(alphaR);
        const T delta2 = square(delta);
        const T delta3 = delta2 * delta;
        const T delta4 = square(delta2);
        return (delta2 - delta3 + delta4 * T(7.0 / 12)) * T(Dr);
    }

    template<Scalar T, class EwaldType>
    T Q_TIP4P<T, EwaldType>::modifiedMorseForce(T r) {
        const T delta = (r - T(equalR)) * T(alphaR);
        const T delta2 = square(delta);
        const T delta3 = delta2 * delta;
        return -(delta * T(2) - delta2 * T(3) + delta3 * T(7.0 / 3)) * T(Dr * alphaR);
    }

    template<Scalar T, class EwaldType>
    auto Q_TIP4P<T, EwaldType>::makeCellWithoutH(const MDCellType& original) -> MDCellType {
        const size_t numMolecule = original.getNumParticle() / 3;
        return MDCellType(original.getLattice(), original.getPos().bottomRows(2 * numMolecule), original.getMassVec());
    }

    template<Scalar T, class EwaldType>
    bool Q_TIP4P<T, EwaldType>::isCellOrdered(const MDCellType& cell) {
        return Base::isCellOrdered(cell, 1, 8);
    }
}

namespace Physica {
    template<Scalar T, class EwaldType>
    class Traits<Q_TIP4P<T, EwaldType>> {
    public:
        constexpr static bool IsPeriodBoundary = true;
        constexpr static bool IsLatticeDependent = true;
        constexpr static bool IsContractable = true;
    };
}
