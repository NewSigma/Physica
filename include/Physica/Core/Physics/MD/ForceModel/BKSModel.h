/*
 * Copyright 2023 WeiBo He.
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

#include "Ewald/RSpaceEwald.h"
#include "ForceModelImpl/AABModel.h"
#include "PairModel.h"

namespace Physica::Core {
    template<class ScalarType, class EwaldType, bool IsSmallCell> class BKSModel;

    namespace Internal {
        template<class T, class EwaldType, bool B>
        class Traits<BKSModel<T, EwaldType, B>> {
        public:
            using ScalarType = T;
            constexpr static bool IsPeriodBoundary = true;
            constexpr static bool IsLatticeDependent = true;
            constexpr static bool IsPotDependOnAtomIndex = true;
            constexpr static bool IsSmallCell = B;
        };
    }
    /**
     * BKS potential for SiO2 or Al-P-O, as introduced in [1]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 64, 1955 (1990); https://doi.org/10.1103/PhysRevLett.64.1955
     */
    template<class ScalarType, class EwaldType, bool IsSmallCell>
    class BKSModel : public PairModel<BKSModel<ScalarType, EwaldType, IsSmallCell>> {
        using This = BKSModel<ScalarType, EwaldType, IsSmallCell>;
        using Base = PairModel<This>;
        using AABModelType = AABModel<ScalarType>;
    public:
        using typename Base::PlainScalar;
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Vector3D;

        constexpr static double angstormInBohr = PhyConst<AU>::angstormToBohr(1);
        constexpr static double angstormInBohr2 = angstormInBohr * angstormInBohr;
        constexpr static double angstormInBohr4 = angstormInBohr2 * angstormInBohr2;
        constexpr static double angstormInBohr6 = angstormInBohr4 * angstormInBohr2;
        constexpr static double A_OO = PhyConst<AU>::eVToHartree(1388.773);
        constexpr static double b_OO = 2.76 * PhyConst<AU>::bohrToAngstorm(1);
        constexpr static double c_OO = PhyConst<AU>::eVToHartree(175) * angstormInBohr6;
        constexpr static double chargeO = -1.2;

        constexpr static double A_SiO = PhyConst<AU>::eVToHartree(18003.7572);
        constexpr static double b_SiO = 4.87318 * PhyConst<AU>::bohrToAngstorm(1);
        constexpr static double c_SiO = PhyConst<AU>::eVToHartree(133.5381) * angstormInBohr6;
        constexpr static double chargeSi = -2 * chargeO;
    private:
        EwaldType ewald;
    public:
        BKSModel(const MDCellType& refer_cell, PlainScalar cutoff, EwaldType ewald_);
        BKSModel(const BKSModel&) = default;
        BKSModel(BKSModel&&) noexcept = default;
        ~BKSModel() = default;
        /* Operators */
        BKSModel& operator=(BKSModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;

        [[nodiscard]] inline ScalarType potentialEnergy(const MDCellType& cell) const;

        template<class Executor> [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const;
        template<class Executor> [[nodiscard]] inline Vector<ScalarType> force_long(const MDCellType& cell) const;

        [[nodiscard]] inline LatticeMatrix virial(const MDCellType& cell) const;
        void swap(BKSModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumMolecule() const noexcept { return getNumParticle() / 3; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return ewald.getNumParticle(); }
        [[nodiscard]] const EwaldType& getEwald() const noexcept { return ewald; }
        [[nodiscard]] const typename MDCellType::LatticeMatrix& getLattice() const noexcept { return ewald.getLattice(); }
        /* Setters */
        inline void setLattice(LatticeMatrix lattice);
        /* Static members */
        inline static PermutationMatrix<ScalarType> sortPosition(MDCellType& cell);
    private:
        [[nodiscard]] inline static bool isCellOrdered(const MDCellType& cell);
    };

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    BKSModel<ScalarType, EwaldType, IsSmallCell>::BKSModel(const MDCellType& refer_cell, PlainScalar cutoff, EwaldType ewald_)
            : Base(cutoff)
            , ewald(std::move(ewald_)) {
        assert(refer_cell.getNumParticle() % 3 == 0 && "[Error]: This is not a cell of SiO2");
        const size_t numMolecule = refer_cell.getNumParticle() / 3;
        ewald = RSpaceEwald<ScalarType, IsSmallCell>(refer_cell.getLattice(), AABModelType::makeCharges(numMolecule, chargeO, chargeSi));
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    ScalarType BKSModel<ScalarType, EwaldType, IsSmallCell>::pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        const size_t numMolecule = getNumMolecule();
        if (i > j)
            std::swap(i, j);
        const bool isSiSi = i >= 2 * numMolecule;
        if (isSiSi)
            return ScalarType(0);
        const bool isOO = j < 2 * numMolecule;
        const double A = isOO ? A_OO : A_SiO;
        const double b = isOO ? b_OO : b_SiO;
        const double c = isOO ? c_OO : c_SiO;
        const ScalarType r4 = square(r2);
        return ScalarType(A) * exp(ScalarType(-b) * r) - ScalarType(c) / (r2 * r4);
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    ScalarType BKSModel<ScalarType, EwaldType, IsSmallCell>::force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        const size_t numMolecule = getNumMolecule();
        if (i > j)
            std::swap(i, j);
        const bool isSiSi = i >= 2 * numMolecule;
        if (isSiSi)
            return ScalarType(0);
        const bool isOO = j < 2 * numMolecule;
        const double A = isOO ? A_OO : A_SiO;
        const double b = isOO ? b_OO : b_SiO;
        const double c = isOO ? c_OO : c_SiO;
        const ScalarType r4 = square(r2);
        return ScalarType(A * b) * exp(ScalarType(-b) * r) - ScalarType(6 * c) / (r * r2 * r4);
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    inline ScalarType BKSModel<ScalarType, EwaldType, IsSmallCell>::potentialEnergy(const MDCellType& cell) const {
        return Base::potentialEnergy(cell) + ewald.potentialEnergy(cell.getPos());
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    template<class Executor>
    Vector<ScalarType> BKSModel<ScalarType, EwaldType, IsSmallCell>::force(const MDCellType& cell) const {
        Vector<ScalarType> result;
        forceAsync<Vector<ScalarType>, Executor>(cell, result);
        return result;
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    template<class VectorType, class Executor>
    void BKSModel<ScalarType, EwaldType, IsSmallCell>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        static_assert(!Internal::Traits<Executor>::isCudaEnabled, "[Error]: Cuda is not supported");
        assert(cell.getNumParticle() % 3 == 0);
        auto future = Executor::schedule([this, &cell, &result]() {
            result = force_short<Executor>(cell);
        });

        const Vector<ScalarType> coulomb = force_long<Executor>(cell);
        Executor::auto_wait(future);
        result += coulomb;
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    template<class Executor>
    Vector<ScalarType> BKSModel<ScalarType, EwaldType, IsSmallCell>::force_short(const MDCellType& cell) const {
        return Base::template force_short<Executor>(cell);
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    template<class Executor>
    inline Vector<ScalarType> BKSModel<ScalarType, EwaldType, IsSmallCell>::force_long(const MDCellType& cell) const {
        return ewald.template force<Executor>(cell.getPos());
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    typename BKSModel<ScalarType, EwaldType, IsSmallCell>::LatticeMatrix
    inline BKSModel<ScalarType, EwaldType, IsSmallCell>::virial(const MDCellType& cell) const {
        return Base::virial(cell) + ewald.virial(cell.getPos());
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    void BKSModel<ScalarType, EwaldType, IsSmallCell>::swap(BKSModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        ewald.swap(obj.ewald);
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    inline void BKSModel<ScalarType, EwaldType, IsSmallCell>::setLattice(LatticeMatrix lattice) {
        ewald.setLattice(std::move(lattice));
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    inline PermutationMatrix<ScalarType> BKSModel<ScalarType, EwaldType, IsSmallCell>::sortPosition(MDCellType& cell) {
        return AABModelType::sortPosition(cell, 8, 14);
    }

    template<class ScalarType, class EwaldType, bool IsSmallCell>
    inline bool BKSModel<ScalarType, EwaldType, IsSmallCell>::isCellOrdered(const MDCellType& cell) {
        return AABModelType::isCellOrdered(cell, 8, 14);
    }
}
