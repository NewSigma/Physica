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

#include <algorithm>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Physics/MD/MDImpl/CellList.h"
#include "EmptyForceModel.h"

namespace Physica::Core {
    /**
     * Member variable pot_shift is referenced from [1]
     * 
     * References:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:205
     */
    template<class Derived>
    class PairModel : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
        using TraitType = Internal::Traits<Derived>;

        constexpr static bool IsPotDependOnAtomIndex = TraitType::IsPotDependOnAtomIndex;
    public:
        using ScalarType = typename TraitType::ScalarType;
        using PosScalarType = typename TraitType::PosScalarType;
        constexpr static unsigned int Dim = 3;

        using PlainScalar = typename ScalarType::PlainScalar;        
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using CellListType = CellList<ScalarType, PosScalarType>;
        using Index3D = typename GridBase::Index3D;
        using Vector3D = Vector<PosScalarType, 3>;
        using ForceConstMatrix = typename EmptyForceModel<ScalarType, PosScalarType, Dim>::ForceConstMatrix;
    private:
        PlainScalar cutoff;
        ScalarType squared_cutoff;
        ScalarType pot_shift;
    public:
        ~PairModel() = default;
        /* Operations */
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType forceConst_functor(ScalarType r, ScalarType r2) const;

        [[nodiscard]] ScalarType potentialEnergy(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const;
        [[nodiscard]] inline ScalarType potentialEnergy(const MDCellType& cell) const;

        template<class Executor, bool IsSmallCell = false>
        [[nodiscard]] Vector<ScalarType> force(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const;
        template<class Executor, bool IsSmallCell = false>
        [[nodiscard]] inline Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor, bool IsSmallCell = false>
        void forceAsync(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos, ContinuousVector<VectorType>& result) const;
        template<class VectorType, class Executor, bool IsSmallCell = false>
        inline void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getNumParticle() * 3, 0); }

        [[nodiscard]] ScalarType forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(PairModel& pair) noexcept;
        /* Getters */
        [[nodiscard]] const PlainScalar& getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] const ScalarType& getSquaredCutoff() const noexcept { return squared_cutoff; }
        /* Setters */
        void setCutoff(PlainScalar cutoff_);
    protected:
        PairModel() = default;
        PairModel(PlainScalar cutoff_);
        PairModel(const PairModel&) = default;
        PairModel(PairModel&&) noexcept = default;
        /* Operators */
        PairModel& operator=(PairModel pair) noexcept;
    private:
        [[nodiscard]] ScalarType forceConstImpl(const MDCellType& cell, size_t dof1, size_t dof2) const;
    };

    template<class Derived>
    PairModel<Derived>::PairModel(PlainScalar cutoff_) {
        setCutoff(std::move(cutoff_));
    }

    template<class Derived>
    PairModel<Derived>& PairModel<Derived>::operator=(PairModel pair) noexcept {
        swap(pair);
        return *this;
    }

    template<class Derived>
    inline typename PairModel<Derived>::ScalarType PairModel<Derived>::pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        return Base::getDerived().pot_functor(i, j, r, r2);
    }

    template<class Derived>
    inline typename PairModel<Derived>::ScalarType PairModel<Derived>::force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        return Base::getDerived().force_functor(i, j, r, r2);
    }

    template<class Derived>
    inline typename PairModel<Derived>::ScalarType PairModel<Derived>::forceConst_functor(ScalarType r, ScalarType r2) const {
        return Base::getDerived().forceConst_functor(r, r2);
    }

    template<class Derived>
    typename PairModel<Derived>::ScalarType PairModel<Derived>::potentialEnergy(
            const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const {
        const auto& pos = cartesianPos;
        const auto range = MDCellType::estimateRange(lattice, cutoff);
        const size_t numParticle = pos.getRow();

        ScalarType result = 0;
        MDCellType::forCellInRange(range, lattice,
            [this, numParticle, &pos, &result](Vector3D delta) {
                ScalarType temp = 0;
                for (size_t i = 0; i < numParticle; ++i) {
                    Vector3D from = pos.row(i) + delta;
                    for (size_t j = i; j < numParticle; ++j) {
                        const ScalarType r2 = (from - pos.row(j)).squaredNorm();
                        const bool isNotSelf = PlainScalar(std::numeric_limits<ScalarType>::min()) < r2;
                        if (isNotSelf && r2 < squared_cutoff) {
                            const ScalarType dist = sqrt(r2);
                            if constexpr (IsPotDependOnAtomIndex)
                                temp += pot_functor(i, j, dist, r2);
                            else
                                temp += pot_functor(i, j, dist, r2) - pot_shift;
                        }
                    }
                }
                result += temp;
            });
        return result;
    }

    template<class Derived>
    inline typename PairModel<Derived>::ScalarType PairModel<Derived>::potentialEnergy(const MDCellType& cell) const {
        return potentialEnergy(cell.getLattice(), cell.getPos());
    }

    template<class Derived>
    template<class Executor, bool IsSmallCell>
    Vector<typename PairModel<Derived>::ScalarType> PairModel<Derived>::force(
            const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const {
        const size_t DOF = cartesianPos.getRow() * cartesianPos.getColumn();
        Vector<ScalarType> result(DOF);
        forceAsync<Vector<ScalarType>, Executor, IsSmallCell>(lattice, cartesianPos, result);
        return result;
    }

    template<class Derived>
    template<class Executor, bool IsSmallCell>
    inline Vector<typename PairModel<Derived>::ScalarType> PairModel<Derived>::force(const MDCellType& cell) const {
        return force<Executor, IsSmallCell>(cell.getLattice(), cell.getPos());
    }

    template<class Derived>
    template<class VectorType, class Executor, bool IsSmallCell>
    void PairModel<Derived>::forceAsync(
            const LatticeMatrix& lattice, const PositionMatrix& cartesianPos, ContinuousVector<VectorType>& result) const {
        static_assert(!Internal::Traits<Executor>::isCudaEnabled, "[Error]: Cuda is not supported");
        const auto& pos = cartesianPos;
        result = ScalarType(0);
        if constexpr (IsSmallCell) {
            const auto range = MDCellType::estimateRange(lattice, cutoff);
            const size_t numParticle = pos.getRow();
            MDCellType::forCellInRange(range, lattice,
                [this, pos, numParticle, &result](Vector3D delta) {
                    Vector3D r, from;
                    for (size_t i = 0; i < numParticle; ++i) {
                        auto force_i = result.segment(3 * i, 3 * i + 3);
                        from = pos.row(i) + delta;
                        for (size_t j = i; j < numParticle; ++j) {
                            auto force_j = result.segment(3 * j, 3 * j + 3);
                            auto to = pos.row(j);
                            r = to.asVector() - from;
                            const ScalarType r2 = r.squaredNorm();
                            const bool isNotSelf = ScalarType(std::numeric_limits<ScalarType>::min()) < r2;
                            if (isNotSelf && r2 < squared_cutoff) {
                                const ScalarType dist = sqrt(r2);
                                const ScalarType f_norm = force_functor(i, j, dist, r2);
                                r *= f_norm / dist;
                                const Vector3D& f = r;
                                force_i -= f;
                                force_j += f;
                            }
                        }
                    }
                });
        }
        else {
            const CellListType cellList(lattice, pos, cutoff);
            Utils::Array<size_t> arr1{};
            cellList.forCellInList([this, pos, &arr1, &result, &cellList](Index3D center) {
                cellList.forAtomInCell(center, [&arr1](size_t atom) {
                    arr1.append(atom);
                });
                /* Current cell */ {
                    const size_t length = arr1.getLength();
                    for (size_t i = 0; i + 1 < length; ++i) {
                        const size_t atom1 = arr1[i];
                        Vector3D f(3, 0);
                        const auto from = pos.row(atom1);
                        for (size_t j = i + 1; j < length; ++j) {
                            const size_t atom2 = arr1[j];
                            const auto to = pos.row(atom2);
                            Vector3D r = to.asVector() - from;
                            const ScalarType r2 = r.squaredNorm();
                            if (r2 < squared_cutoff) {
                                const ScalarType dist = sqrt(r2);
                                const ScalarType f_norm = force_functor(atom1, atom2, dist, r2);
                                r *= ScalarType(f_norm / dist);
                                f -= r;
                                auto f2 = result.template segment<3>(3 * atom2, 3 * atom2 + 3);
                                f2 += r;
                            }
                        }
                        auto f1 = result.template segment<3>(3 * atom1, 3 * atom1 + 3);
                        f1 += f;
                    }
                }
                Utils::Array<size_t> arr2{};
                cellList.forReducedNeighInRange(center, [this, pos, &arr1, &arr2, &result, &cellList](Vector3D translate, Index3D neigh) {
                    cellList.forAtomInCell(neigh, [&arr2](size_t atom) {
                        arr2.append(atom);
                    });
                    std::sort(arr2.begin(), arr2.end());

                    for (const size_t atom1 : arr1) {
                        const Vector3D from = pos.row(atom1) - translate;
                        Vector3D r, f(3, 0);
                        for (const size_t atom2 : arr2) {
                            const auto to = pos.row(atom2);
                            r = to.asVector() - from;
                            const ScalarType r2 = r.squaredNorm();
                            if (r2 < squared_cutoff) {
                                const ScalarType dist = sqrt(r2);
                                auto f2 = result.template segment<3>(3 * atom2, 3 * atom2 + 3);
                                const ScalarType f_norm = force_functor(atom1, atom2, dist, r2);
                                r *= ScalarType(f_norm / dist);
                                f -= r;
                                f2 += r;
                            }
                        }
                        auto f1 = result.template segment<3>(3 * atom1, 3 * atom1 + 3);
                        f1 += f;
                    }
                    arr2.clear();
                });
                arr1.clear();
            });
        }
    }

    template<class Derived>
    template<class VectorType, class Executor, bool IsSmallCell>
    inline void PairModel<Derived>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        forceAsync<VectorType, Executor, IsSmallCell>(cell.getLattice(), cell.getPos(), result);
    }

    template<class Derived>
    typename PairModel<Derived>::ScalarType PairModel<Derived>::forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const {
        const size_t atom1 = dof1 / 3U;
        const size_t atom2 = dof2 / 3U;
        [[unlikely]] if (atom1 == atom2) {
            const size_t dir2 = dof2 % 3U;
            ScalarType result = 0;
            for (size_t i = 0; i < cell.getNumParticle(); ++i) {
                if (i == atom2)
                    continue;
                const size_t dof_i = 3 * i + dir2;
                result += forceConstImpl(cell, dof1, dof_i);
            }
            return result;
        }
        return -forceConstImpl(cell, dof1, dof2);
    }

    template<class Derived>
    typename PairModel<Derived>::ForceConstMatrix PairModel<Derived>::forceConst(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const size_t dof = cell.getDOF();
        ForceConstMatrix result(dof);
        for (size_t atom_r = 0; atom_r < numParticle; ++atom_r) {
            for (size_t dim_r = 0; dim_r < 3; ++dim_r) {
                const size_t r = atom_r * 3 + dim_r;
                for (size_t c = (atom_r + 1) * Dim; c < dof; ++c)
                    result(r, c) = forceConst(cell, r, c);

                //Handle self interaction
                for (size_t dim_c = dim_r; dim_c < 3; ++dim_c) {
                    const size_t c = atom_r * 3 + dim_c;
                    ScalarType temp = 0;
                    for (size_t atom_c = 0; atom_c < numParticle; ++atom_c) {
                        if (atom_r == atom_c)
                            continue;
                        temp += result(r, atom_c * 3 + dim_c);
                    }
                    result(r, c) = -temp;
                }
            }
        }
        return result;
    }
    /**
     * Reference:
     * [1] M. J. Louwerse and E. J. Baerends, Chem. Phys. Lett. 421, 138 (2006).
     */
    template<class Derived>
    typename PairModel<Derived>::LatticeMatrix
    PairModel<Derived>::virial(const MDCellType& cell) const {
        const auto& pos = cell.getPos();
        const size_t numParticle = cell.getNumParticle();
        const CellListType cellList(cell, cutoff);

        LatticeMatrix result(3, 3, 0);
        for (size_t atom1 = 0; atom1 < numParticle; ++atom1) {
            const Index3D center = cellList.getAtomCellMap()[atom1];
            cellList.forNeighInRange(center, [this, atom1, pos, &cellList, &result](Vector3D translate, Index3D neigh) {
                const Vector3D from = pos.row(atom1) - translate;
                const auto& subCell = cellList(neigh);

                Vector3D r, f;
                auto ite = subCell.cbegin();
                bool isValid = ite != subCell.cend();
                size_t atomToSolve;
                if (isValid)
                    atomToSolve = *(ite++);

                while (isValid) {
                    const size_t atom2 = atomToSolve;
                    isValid = ite != subCell.cend();
                    if (isValid)
                        atomToSolve = *(ite++);

                    const bool isSelf = atom1 == atom2;
                    if (isSelf)
                        continue;
                    const auto to = pos.row(atom2);
                    r = to.asVector() - from;
                    const ScalarType r2 = r.squaredNorm();
                    if (r2 < squared_cutoff) {
                        const ScalarType dist = sqrt(r2);
                        const ScalarType f_norm = force_functor(atom1, atom2, dist, r2);
                        f = r * ScalarType(-f_norm / dist);
                        result += f * r.transpose();
                    }
                }
            });
        }
        result *= reciprocal(ScalarType(cell.getVolume() * 2.0));
        return result;
    }

    template<class Derived>
    void PairModel<Derived>::swap(PairModel& pair) noexcept {
        assert(this != &pair && "[Error]: Self swap is likely a bug");
        cutoff.swap(pair.cutoff);
        squared_cutoff.swap(pair.squared_cutoff);
        pot_shift.swap(pair.pot_shift);
    }

    template<class Derived>
    void PairModel<Derived>::setCutoff(PlainScalar cutoff_) {
        cutoff = std::move(cutoff_);
        squared_cutoff = square(cutoff);
        if constexpr (!IsPotDependOnAtomIndex)
            pot_shift = pot_functor(0, 0, cutoff, squared_cutoff);
    }

    template<class Derived>
    typename PairModel<Derived>::ScalarType PairModel<Derived>::forceConstImpl(const MDCellType& cell, size_t dof1, size_t dof2) const {
        static_assert(!IsPotDependOnAtomIndex, "[Error]: It is assumed force not depends on atom index");
        constexpr int unused = 0;
        const size_t atom1 = dof1 / 3U;
        const size_t atom2 = dof2 / 3U;
        assert(atom1 != atom2 && "[Error]: The function is not responsible for this case");
        const size_t dir1 = dof1 % 3U;
        const size_t dir2 = dof2 % 3U;
        const Vector3D delta = cell.getPos().row(atom2).asVector() - cell.getPos().row(atom1).asVector();
        const ScalarType squaredNorm = delta.squaredNorm();
        const ScalarType norm = sqrt(squaredNorm);
        const ScalarType factor = delta[dir1] * delta[dir2] / squaredNorm;
        const ScalarType term1 = forceConst_functor(norm, squaredNorm) * factor;
        const ScalarType term2 = force_functor(unused, unused, norm, squaredNorm) / norm * (ScalarType(dir1 == dir2 ? 1.0 : 0.0) - factor);
        return term1 - term2;
    }
}
