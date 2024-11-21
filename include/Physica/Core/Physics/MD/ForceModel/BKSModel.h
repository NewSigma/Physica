/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/PlainStruct.h"
#include "Ewald/RSpaceEwald.h"
#include "ForceModelImpl/AABModel.h"
#include "PairModel.h"

namespace Physica::Core {
    /**
     * BKS potential for SiO2 or Al-P-O, as introduced in [1]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 64, 1955 (1990); https://doi.org/10.1103/PhysRevLett.64.1955
     */
    template<Scalar T, class EwaldType, bool AvoidTooNear>
    class BKSModel : public PairModel<BKSModel<T, EwaldType, AvoidTooNear>> {
        using This = BKSModel<T, EwaldType, AvoidTooNear>;
        using Base = PairModel<This>;
        using AABModelType = AABModel<T>;
        using REwaldType = Traits<EwaldType>::REwaldType;
        using DoublePotType = std::conditional<AvoidTooNear, T, PlainStruct<void>>::type;
    public:
        using typename Base::ValueType;
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;

        constexpr static double angstormInBohr = PhyConst<AU>::angstormToBohr(1);
        constexpr static double angstormInBohr2 = angstormInBohr * angstormInBohr;
        constexpr static double angstormInBohr4 = angstormInBohr2 * angstormInBohr2;
        constexpr static double angstormInBohr6 = angstormInBohr4 * angstormInBohr2;
        constexpr static double A_OO = PhyConst<AU>::eVToHartree(1388.773);
        constexpr static double b_OO = 2.76 * PhyConst<AU>::bohrToAngstorm(1);
        constexpr static double c_OO = PhyConst<AU>::eVToHartree(175) * angstormInBohr6;
        constexpr static double r0_OO = 2.8414313142730038;
        constexpr static double chargeO = -1.2;

        constexpr static double A_SiO = PhyConst<AU>::eVToHartree(18003.7572);
        constexpr static double b_SiO = 4.87318 * PhyConst<AU>::bohrToAngstorm(1);
        constexpr static double c_SiO = PhyConst<AU>::eVToHartree(133.5381) * angstormInBohr6;
        constexpr static double r0_SiO = 2.075401485810223;
        constexpr static double chargeSi = -2 * chargeO;
    private:
        EwaldType ewald;
        [[no_unique_address]] DoublePotType doublePot0_OO;
        [[no_unique_address]] DoublePotType doublePot0_SiO;
    public:
        BKSModel(const MDCellType& refer_cell, ValueType cutoff, EwaldType ewald_);
        BKSModel(const BKSModel&) = default;
        BKSModel(BKSModel&&) noexcept = default;
        ~BKSModel() = default;
        /* Operators */
        BKSModel& operator=(BKSModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] T force_functor(size_t i, size_t j, T r, T r2) const;

        [[nodiscard]] inline T potentialV(const MDCellType& cell) const;

        template<class Executor> [[nodiscard]] VectorND<T> force(const MDCellType& cell);
        template<Vector V, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<V>& result);
        template<class Executor> [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const;
        template<class Executor> [[nodiscard]] inline VectorND<T> force_long(const MDCellType& cell);

        [[nodiscard]] inline LatticeMatrix virial(const MDCellType& cell);
        void swap(BKSModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return getNumParticle() / 3; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return ewald.getNumParticle(); }
        [[nodiscard]] const EwaldType& getEwald() const noexcept { return ewald; }
        [[nodiscard]] const MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Setters */
        inline void setLattice(LatticeMatrix lattice);
        /* Static members */
        inline static PermutationMatrix<T> sortPosition(MDCellType& cell);
    private:
        [[nodiscard]] inline static T pot_functor_impl(T A, T b, T c, T r, T r2);
        [[nodiscard]] inline static bool isCellOrdered(const MDCellType& cell);
    };

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    BKSModel<T, EwaldType, AvoidTooNear>::BKSModel(const MDCellType& refer_cell, ValueType cutoff, EwaldType ewald_)
            : Base(cutoff)
            , ewald(std::move(ewald_)) {
        assert(refer_cell.getNumParticle() % 3 == 0 && "[Error]: This is not a cell of SiO2");
        const size_t numMolecule = refer_cell.getNumParticle() / 3;
        ewald = REwaldType(refer_cell.getLattice(), AABModelType::makeCharges(numMolecule, chargeO, chargeSi));
        if constexpr (AvoidTooNear) {
            doublePot0_OO = pot_functor_impl(A_OO, b_OO, c_OO, r0_OO, r0_OO * r0_OO);
            doublePot0_SiO = pot_functor_impl(A_SiO, b_SiO, c_SiO, r0_SiO, r0_SiO * r0_SiO);
        }
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    T BKSModel<T, EwaldType, AvoidTooNear>::pot_functor(size_t i, size_t j, T r, T r2) const {
        const size_t numMolecule = getNumMolecule();
        if (i > j)
            std::swap(i, j);
        const bool isSiSi = i >= 2 * numMolecule;
        if (isSiSi)
            return T(0);
        const bool isOO = j < 2 * numMolecule;
        const double A = isOO ? A_OO : A_SiO;
        const double b = isOO ? b_OO : b_SiO;
        const double c = isOO ? c_OO : c_SiO;
        const T pot = pot_functor_impl(A, b, c, r, r2);
        if constexpr (AvoidTooNear) {
            const double r0 = isOO ? r0_OO : r0_SiO;
            const T doublePot0 = isOO ? doublePot0_OO : doublePot0_SiO;
            return r < r0 ? (doublePot0 - pot) : pot;
        }
        else
            return pot;
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    T BKSModel<T, EwaldType, AvoidTooNear>::force_functor(size_t i, size_t j, T r, T r2) const {
        const size_t numMolecule = getNumMolecule();
        if (i > j)
            std::swap(i, j);
        const bool isSiSi = i >= 2 * numMolecule;
        if (isSiSi)
            return T(0);
        const bool isOO = j < 2 * numMolecule;
        const double A = isOO ? A_OO : A_SiO;
        const double b = isOO ? b_OO : b_SiO;
        const double c = isOO ? c_OO : c_SiO;
        const T r4 = square(r2);
        const T f = T(A * b) * exp(T(-b) * r) - T(6 * c) / (r * r2 * r4);
        if constexpr (AvoidTooNear) {
            const double r0 = isOO ? r0_OO : r0_SiO;
            return r < r0 ? -f : f;
        }
        else
            return f;
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    inline T BKSModel<T, EwaldType, AvoidTooNear>::potentialV(const MDCellType& cell) const {
        return Base::potentialV(cell) + ewald.potentialV(cell.getPos());
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    template<class Executor>
    VectorND<T> BKSModel<T, EwaldType, AvoidTooNear>::force(const MDCellType& cell) {
        VectorND<T> result;
        forceAsync<VectorND<T>, Executor>(cell, result);
        return result;
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    template<Vector V, class Executor>
    void BKSModel<T, EwaldType, AvoidTooNear>::forceAsync(const MDCellType& cell, ContinuousVector<V>& result) {
        static_assert(!Traits<Executor>::UseCUDA, "[Error]: CUDA is not supported");
        assert(cell.getNumParticle() % 3 == 0);
        auto future = Executor::schedule([this, &cell, &result]() {
            result = force_short<Executor>(cell);
        });

        const VectorND<T> coulomb = force_long<Executor>(cell);
        Executor::auto_wait(future);
        result += coulomb;
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    template<class Executor>
    VectorND<T> BKSModel<T, EwaldType, AvoidTooNear>::force_short(const MDCellType& cell) const {
        return Base::template force_short<Executor>(cell);
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    template<class Executor>
    inline VectorND<T> BKSModel<T, EwaldType, AvoidTooNear>::force_long(const MDCellType& cell) {
        return ewald.template force<Executor>(cell.getPos());
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    BKSModel<T, EwaldType, AvoidTooNear>::LatticeMatrix
    inline BKSModel<T, EwaldType, AvoidTooNear>::virial(const MDCellType& cell) {
        return Base::virial(cell) + ewald.virial(cell.getPos());
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    void BKSModel<T, EwaldType, AvoidTooNear>::swap(BKSModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        ewald.swap(obj.ewald);
        if constexpr (AvoidTooNear) {
            doublePot0_OO.swap(obj.doublePot0_OO);
            doublePot0_SiO.swap(obj.doublePot0_SiO);
        }
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    inline void BKSModel<T, EwaldType, AvoidTooNear>::setLattice(LatticeMatrix lattice) {
        ewald.setLattice(std::move(lattice));
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    inline PermutationMatrix<T> BKSModel<T, EwaldType, AvoidTooNear>::sortPosition(MDCellType& cell) {
        return AABModelType::sortPosition(cell, 8, 14);
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    inline T BKSModel<T, EwaldType, AvoidTooNear>::pot_functor_impl(
            T A, T b, T c, T r, T r2) {
        const T r4 = square(r2);
        return A * exp(-b * r) - c / (r2 * r4);
    }

    template<Scalar T, class EwaldType, bool AvoidTooNear>
    inline bool BKSModel<T, EwaldType, AvoidTooNear>::isCellOrdered(const MDCellType& cell) {
        return AABModelType::isCellOrdered(cell, 8, 14);
    }
}

namespace Physica {
    template<Scalar T, class EwaldType, bool AvoidTooNear>
    class Traits<Core::BKSModel<T, EwaldType, AvoidTooNear>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPeriodBoundary = true;
        constexpr static bool IsLatticeDependent = true;
        constexpr static bool IsPotDependOnAtomIndex = true;
        constexpr static bool IsSmallCell = Traits<EwaldType>::IsSmallCell;
        constexpr static bool IsContractable = false;
    };
}
