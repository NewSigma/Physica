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
#include "Physica/Core/Physics/MD/MDImpl/CellList.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:205
     */
    template<class ScalarType, class PosScalarType, class PairFunctor>
    class PairModel final {
        using ResultType = typename std::invoke_result<PairFunctor, PosScalarType>::type;
        static_assert(is_scalar<ScalarType>::value && is_scalar<PosScalarType>::value, "[Error]: Invalid ScalarType");
        static_assert(std::is_same<ScalarType, ResultType>::value, "[Error]: Invalid PairFunctor");
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using CellListType = CellList<ScalarType, PosScalarType>;
        using Index3D = typename CellListType::Index3D;
        using Vector3D = Vector<PosScalarType, 3>;
    private:
        ScalarType cutoff;
        ScalarType squared_cutoff;
        ScalarType pot_shift;
        PairFunctor force_functor;
        PairFunctor pot_functor;
    public:
        PairModel(ScalarType cutoff_, PairFunctor functor_, PairFunctor pot_functor_);
        PairModel(const PairModel&) = default;
        PairModel(PairModel&&) noexcept = default;
        ~PairModel() = default;
        /* Operators */
        PairModel& operator=(PairModel pair) noexcept;
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getNumParticle() * 3, 0); }
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(PairModel& pair) noexcept;
    };

    template<class ScalarType, class PosScalarType, class PairFunctor>
    PairModel<ScalarType, PosScalarType, PairFunctor>::PairModel(ScalarType cutoff_, PairFunctor functor_, PairFunctor pot_functor_)
            : cutoff(std::move(cutoff_))
            , force_functor(std::move(functor_))
            , pot_functor(std::move(pot_functor_)) {
        squared_cutoff = square(cutoff);
        pot_shift = pot_functor(cutoff);
    }

    template<class ScalarType, class PosScalarType, class PairFunctor>
    PairModel<ScalarType, PosScalarType, PairFunctor>& PairModel<ScalarType, PosScalarType, PairFunctor>::operator=(PairModel pair) noexcept {
        swap(pair);
        return *this;
    }

    template<class ScalarType, class PosScalarType, class PairFunctor>
    template<class Executor>
    Vector<ScalarType> PairModel<ScalarType, PosScalarType, PairFunctor>::force(const MDCellType& cell) const {
        const auto& pos = cell.getPos();
        const size_t numParticle = cell.getNumParticle();
        const CellListType cellList(cell, cutoff);

        Vector<ScalarType> force(3 * cell.getNumParticle(), 0);
        for (size_t atom1 = 0; atom1 < numParticle; ++atom1) {
            const Index3D center = cellList.getAtomCellMap()[atom1];
            cellList.forNeighInRange(center, [this, atom1, pos, &cellList, &force](Vector3D translate, Index3D neigh) {
                const Vector3D from = pos.row(atom1) - translate;
                const auto& subCell = cellList(neigh);

                Vector3D r, f(3, 0);
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

                    const bool isDoubleCounted = atom2 < atom1;
                    if (isDoubleCounted)
                        continue;
                    const auto to = pos.row(atom2);
                    auto f2 = force.segment(3 * atom2, 3 * atom2 + 3);
                    r = to.asVector() - from;
                    const ScalarType r2 = r.squaredNorm();
                    const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                    if (isNotSelf && r2 < squared_cutoff) {
                        const ScalarType dist = sqrt(r2);
                        const ScalarType f_norm = force_functor(dist);
                        r *= ScalarType(f_norm / dist);
                        f -= r;
                        f2 += r;
                    }
                }
                auto f1 = force.segment(3 * atom1, 3 * atom1 + 3);
                f1 += f;
            });
        }
        return force;
    }

    template<class ScalarType, class PosScalarType, class PairFunctor>
    ScalarType PairModel<ScalarType, PosScalarType, PairFunctor>::potentialEnergy(const MDCellType& cell) const {
        const auto& pos = cell.getPos();
        const auto range = MDCellType::estimateRange(cell.getLattice(), cutoff);
        const size_t numParticle = cell.getNumParticle();

        ScalarType result = 0;
        MDCellType::forCellInRange(range, cell.getLattice(),
            [this, numParticle, &pos, &result](Vector3D delta) {
                ScalarType temp = 0;
                for (size_t i = 0; i < numParticle; ++i) {
                    Vector3D from = pos.row(i) + delta;
                    for (size_t j = i; j < numParticle; ++j) {
                        const ScalarType r2 = (from - pos.row(j)).squaredNorm();
                        const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                        if (isNotSelf && r2 < squared_cutoff) {
                            const ScalarType dist = sqrt(r2);
                            temp += pot_functor(dist) - pot_shift;
                        }
                    }
                }
                result += temp;
            });
        return result;
    }

    template<class ScalarType, class PosScalarType, class PairFunctor>
    typename PairModel<ScalarType, PosScalarType, PairFunctor>::LatticeMatrix
    PairModel<ScalarType, PosScalarType, PairFunctor>::virial(const MDCellType& cell) const {
        const auto& pos = cell.getPos();
        const size_t numParticle = cell.getNumParticle();
        const CellListType cellList(cell, cutoff);

        LatticeMatrix result(3, 3, 0);
        for (size_t atom1 = 0; atom1 < numParticle; ++atom1) {
            const Index3D center = cellList.getAtomCellMap()[atom1];
            Vector3D f1(3, 0);
            cellList.forNeighInRange(center, [this, atom1, pos, &cellList, &f1](Vector3D translate, Index3D neigh) {
                const Vector3D from = pos.row(atom1) - translate;
                const auto& subCell = cellList(neigh);

                Vector3D r, f(3, 0);
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
                        const ScalarType f_norm = force_functor(dist);
                        r *= ScalarType(f_norm / dist);
                        f -= r;
                    }
                }
                f1 += f;
            });
            result += f1 * pos.row(atom1).asVector().transpose();
        }
        result *= reciprocal(ScalarType(cell.getVolume() * 2.0));
        return result;
    }

    template<class ScalarType, class PosScalarType, class PairFunctor>
    void PairModel<ScalarType, PosScalarType, PairFunctor>::swap(PairModel& pair) noexcept {
        cutoff.swap(pair.cutoff);
        squared_cutoff.swap(squared_cutoff);
        pot_shift.swap(pot_shift);
        std::swap(force_functor, force_functor);
        std::swap(pot_functor, pot_functor);
    }
}
