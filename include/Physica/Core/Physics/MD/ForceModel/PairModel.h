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

namespace Physica::Core {
    /**
     * References:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:205
     */
    template<class Derived>
    class PairModel : public Utils::CRTPBase<Derived> {
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        using PosScalarType = typename Internal::Traits<Derived>::PosScalarType;
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using CellListType = CellList<ScalarType, PosScalarType>;
        using Index3D = typename CellListType::Index3D;
        using Vector3D = Vector<PosScalarType, 3>;
        
    private:
        static_assert(is_scalar<ScalarType>::value && is_scalar<PosScalarType>::value, "[Error]: Invalid ScalarType");
        using Base = Utils::CRTPBase<Derived>;

        ScalarType cutoff;
        ScalarType squared_cutoff;
        ScalarType pot_shift;
    public:
        PairModel(ScalarType cutoff_);
        PairModel(const PairModel&) = default;
        PairModel(PairModel&&) noexcept = default;
        ~PairModel() = default;
        /* Operators */
        PairModel& operator=(PairModel pair) noexcept;
        /* Operations */
        [[nodiscard]] ScalarType force_functor(ScalarType r, ScalarType r2) const { return Base::getDerived().force_functor(r, r2); }
        [[nodiscard]] ScalarType pot_functor(ScalarType r, ScalarType r2) const { return Base::getDerived().pot_functor(r, r2); }
        template<class Executor, bool IsSmallCell = false> [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor> [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getNumParticle() * 3, 0); }
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(PairModel& pair) noexcept;
    };

    template<class Derived>
    PairModel<Derived>::PairModel(ScalarType cutoff_)
            : cutoff(std::move(cutoff_)) {
        squared_cutoff = square(cutoff);
        pot_shift = pot_functor(cutoff, squared_cutoff);
    }

    template<class Derived>
    PairModel<Derived>& PairModel<Derived>::operator=(PairModel pair) noexcept {
        swap(pair);
        return *this;
    }

    template<class Derived>
    template<class Executor, bool IsSmallCell>
    Vector<typename PairModel<Derived>::ScalarType> PairModel<Derived>::force(const MDCellType& cell) const {
        const auto& pos = cell.getPos();
        Vector<ScalarType> force(3 * cell.getNumParticle(), 0);
        if constexpr (IsSmallCell) {
            const auto& lattice = cell.getLattice();
            const auto range = MDCellType::estimateRange(lattice, cutoff);
            const size_t numParticle = cell.getNumParticle();

            MDCellType::forCellInRange(range, lattice,
                [this, pos, numParticle, &force](Vector3D delta) {
                    Vector3D r, from;
                    for (size_t i = 0; i < numParticle; ++i) {
                        auto force_i = force.segment(3 * i, 3 * i + 3);
                        from = pos.row(i) + delta;
                        for (size_t j = i; j < numParticle; ++j) {
                            auto force_j = force.segment(3 * j, 3 * j + 3);
                            auto to = pos.row(j);
                            r = to.asVector() - from;
                            const ScalarType r2 = r.squaredNorm();
                            const bool isNotSelf = ScalarType(std::numeric_limits<ScalarType>::min()) < r2;
                            if (isNotSelf && r2 < squared_cutoff) {
                                const ScalarType dist = sqrt(r2);
                                const ScalarType f_norm = force_functor(dist, r2);
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
            const CellListType cellList(cell, cutoff);
            Utils::Array<size_t> arr1{};
            cellList.forCellInList([this, pos, &arr1, &force, &cellList](Index3D center) {
                for (size_t i : cellList(center))
                    arr1.append(i);
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
                                const ScalarType f_norm = force_functor(dist, r2);
                                r *= ScalarType(f_norm / dist);
                                f -= r;
                                auto f2 = force.template segment<3>(3 * atom2, 3 * atom2 + 3);
                                f2 += r;
                            }
                        }
                        auto f1 = force.template segment<3>(3 * atom1, 3 * atom1 + 3);
                        f1 += f;
                    }
                }
                Utils::Array<size_t> arr2{};
                cellList.forReducedNeighInRange(center, [this, pos, &arr1, &arr2, &force, &cellList](Vector3D translate, Index3D neigh) {
                    for (size_t j : cellList(neigh))
                        arr2.append(j);
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
                                auto f2 = force.template segment<3>(3 * atom2, 3 * atom2 + 3);
                                const ScalarType f_norm = force_functor(dist, r2);
                                r *= ScalarType(f_norm / dist);
                                f -= r;
                                f2 += r;
                            }
                        }
                        auto f1 = force.template segment<3>(3 * atom1, 3 * atom1 + 3);
                        f1 += f;
                    }
                    arr2.clear();
                });
                arr1.clear();
            });
        }
        return force;
    }

    template<class Derived>
    typename PairModel<Derived>::ScalarType PairModel<Derived>::potentialEnergy(const MDCellType& cell) const {
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
                        const bool isNotSelf = ScalarType(std::numeric_limits<ScalarType>::min()) < r2;
                        if (isNotSelf && r2 < squared_cutoff) {
                            const ScalarType dist = sqrt(r2);
                            temp += pot_functor(dist, r2) - pot_shift;
                        }
                    }
                }
                result += temp;
            });
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
                        const ScalarType f_norm = force_functor(dist, r2);
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
        cutoff.swap(pair.cutoff);
        squared_cutoff.swap(squared_cutoff);
        pot_shift.swap(pot_shift);
    }
}
