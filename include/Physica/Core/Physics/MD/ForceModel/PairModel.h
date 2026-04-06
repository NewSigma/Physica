/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Physics/MD/MDImpl/CellList.h"
#include "EmptyForceModel.h"

namespace Physica {
    /**
     * Member variable pot_shift is referenced from [1]
     * 
     * References:
     * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:205
     */
    template<class Derived>
    class PairModel : public CRTPBase<PairModel<Derived>> {
        using This = PairModel<Derived>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;

        constexpr static bool IsPotDependOnAtomIndex = TraitsType::IsPotDependOnAtomIndex;
        constexpr static bool IsSmallCell = TraitsType::IsSmallCell;
    public:
        constexpr static unsigned int Dim = 3;
        using ScalarType = TraitsType::ScalarType;
        using T = ScalarType;
        using Tv = T::ValueType;
        using MDCellType = MDCell<T, Dim>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using CellListType = CellList<Tv>;
        using Vec3D = Vector3D<T>;
        using ForceConstMatrix = EmptyForceModel<T, Dim>::ForceConstMatrix;
    private:
        Tv cutoff;
        Tv squared_cutoff;
        Tv pot_shift;
    public:
        ~PairModel() = default;
        /* Operations */
        [[nodiscard]] CoDiff<T> pot_functor(size_t i, size_t j, const T& r, const T& r2) const;
        [[nodiscard]] CoDiff<T> force_functor(size_t i, size_t j, const T& r, const T& r2) const;
        [[nodiscard]] CoDiff<T> forceConst_functor(const T& r, const T& r2) const;

        [[nodiscard]] CoDiff<T> potentialV(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const;
        [[nodiscard]] CoDiff<T> potentialV(const MDCellType& cell) const;

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<ExecutePolicy P>
        void forceAsync(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos, Vector auto& result) const;
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result) const;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<P>(cell); }
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getNumParticle() * 3, 0); }

        [[nodiscard]] CoDiff<T> forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] CoDiff<ForceConstMatrix> forceConst(const MDCellType& cell) const;

        [[nodiscard]] CoDiff<LatticeMatrix> virial(const MDCellType& cell) const;
        [[nodiscard]] CoDiff<LatticeMatrix> virial(const LatticeMatrix& lattice, const PositionMatrix& pos) const;
        void swap(PairModel& __restrict pair) noexcept;
        /* Getters */
        [[nodiscard]] const Tv& getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] const Tv& getSquaredCutoff() const noexcept { return squared_cutoff; }
        /* Setters */
        void setCutoff(Tv cutoff_);
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return true; }
    protected:
        PairModel() = default;
        PairModel(Tv cutoff_);
        PairModel(const PairModel&) = default;
        PairModel(PairModel&&) noexcept = default;
        /* Operators */
        PairModel& operator=(PairModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void forPairInCutoff(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos, std::invocable<int, int, Vec3D, T, T> auto fn) const;
    private:
        [[nodiscard]] CoDiff<T> forceConstImpl(const MDCellType& cell, size_t dof1, size_t dof2) const;
    };

    template<class Derived>
    PairModel<Derived>::PairModel(Tv cutoff_) {
        setCutoff(std::move(cutoff_));
    }

    template<class Derived>
    auto PairModel<Derived>::pot_functor(size_t i, size_t j, const T& r, const T& r2) const -> CoDiff<T> {
        return Base::getDerived().pot_functor(i, j, r, r2);
    }

    template<class Derived>
    auto PairModel<Derived>::force_functor(size_t i, size_t j, const T& r, const T& r2) const -> CoDiff<T> {
        return Base::getDerived().force_functor(i, j, r, r2);
    }

    template<class Derived>
    auto PairModel<Derived>::forceConst_functor(const T& r, const T& r2) const -> CoDiff<T> {
        return Base::getDerived().forceConst_functor(r, r2);
    }

    template<class Derived>
    auto PairModel<Derived>::potentialV(const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const -> CoDiff<T> {
        const auto& pos = cartesianPos;
        const auto range = MDCellType::estimateRange(lattice, cutoff);
        const size_t numParticle = pos.getRow();

        T result = 0;
        MDCellType::forCellInRange(range, lattice,
            [this, numParticle, &pos, &result](Vec3D delta) {
                T temp = 0;
                for (size_t i = 0; i < numParticle; ++i) {
                    CoDiff<Vec3D> from = pos.row(i) + delta;
                    for (size_t j = i; j < numParticle; ++j) {
                        const T r2 = (from - pos.row(j)).squaredNorm();
                        const bool isNotSelf = !r2.isSubNormal();
                        if (isNotSelf && r2 < squared_cutoff) {
                            const T dist = sqrt(r2);
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
    auto PairModel<Derived>::potentialV(const MDCellType& cell) const -> CoDiff<T> {
        return potentialV(cell.getLattice(), cell.getPos());
    }

    template<class Derived>
    template<ExecutePolicy P>
    VectorND<typename PairModel<Derived>::T> PairModel<Derived>::force(
            const LatticeMatrix& lattice, const PositionMatrix& cartesianPos) const {
        const size_t DOF = cartesianPos.getRow() * cartesianPos.getCol();
        VectorND<T> result(DOF);
        forceAsync<P>(lattice, cartesianPos, result);
        return result;
    }

    template<class Derived>
    template<ExecutePolicy P>
    VectorND<typename PairModel<Derived>::T> PairModel<Derived>::force(const MDCellType& cell) const {
        return force<P>(cell.getLattice(), cell.getPos());
    }

    template<class Derived>
    template<ExecutePolicy P>
    void PairModel<Derived>::forceAsync(
            const LatticeMatrix& lattice, const PositionMatrix& cartesianPos, Vector auto& result) const {
        result = T(0);
        auto kernel = [this, &result](size_t i, size_t j, Vec3D r, const T& norm1, const T& norm2) {
            const T f_norm = force_functor(i, j, norm1, norm2);
            r *= f_norm / norm1;
            auto force_i = result.template segment<Dim>(Dim * i, Dim * i + Dim);
            auto force_j = result.template segment<Dim>(Dim * j, Dim * j + Dim);
            force_i -= r;
            force_j += r;
        };
        forPairInCutoff(lattice, cartesianPos, kernel);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void PairModel<Derived>::forceAsync(const MDCellType& cell, Vector auto& result) const {
        forceAsync<P>(cell.getLattice(), cell.getPos(), result);
    }

    template<class Derived>
    auto PairModel<Derived>::forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const -> CoDiff<T> {
        static_assert(isPeriodBoundary(), "[Error]: Fixed boundary is not implemented");
        const size_t atom1 = dof1 / 3U;
        const size_t atom2 = dof2 / 3U;
        if (atom1 == atom2) [[unlikely]] {
            const size_t dir2 = dof2 % 3U;
            if constexpr (ReverseDiff<T>) {
                Tv result = 0;
                size_t i = 0;
                auto _ = co_for([&]{ return i < cell.getNumParticle(); }, [&]{ ++i; }, [&]() -> CoDiff<T> {
                    if (i == atom2)
                        return {};
                    const size_t dof_i = Dim * i + dir2;
                    const auto fc = forceConstImpl(cell, dof1, dof_i);
                    result += fc;
                    return fc;
                });
                std::ignore = co_yield result;
            }
            else {
                T result = 0;
                for (size_t i = 0; i < cell.getNumParticle(); ++i) {
                    if (i == atom2)
                        continue;
                    const size_t dof_i = Dim * i + dir2;
                    result += forceConstImpl(cell, dof1, dof_i);
                }
                co_return std::move(result);
            }
        }

        const auto result = -forceConstImpl(cell, dof1, dof2);
        if constexpr (ReverseDiff<T>) {
            auto& r = co_yield result.value();
            result.reverse(r.grad());
        }
        else
            co_return std::move(result);
    }

    template<class Derived>
    auto PairModel<Derived>::forceConst(const MDCellType& cell) const -> CoDiff<ForceConstMatrix> {
        using MatrixType = std::conditional<ReverseDiff<T>, DenseSymmMatrix<Tv>, ForceConstMatrix>::type;
        using T1 = MatrixType::ScalarType;
        const size_t numParticle = cell.getNumParticle();
        const size_t dof = cell.getDOF();
        
        MatrixType result(dof);
        for (size_t atom_r = 0; atom_r < numParticle; ++atom_r) {
            for (size_t dim_r = 0; dim_r < Dim; ++dim_r) {
                const size_t r = atom_r * Dim + dim_r;
                for (size_t c = (atom_r + 1) * Dim; c < dof; ++c)
                    result[r, c] = T1(forceConst(cell, r, c));

                //Handle self interaction
                for (size_t dim_c = dim_r; dim_c < Dim; ++dim_c) {
                    const size_t c = atom_r * Dim + dim_c;
                    T temp = 0;
                    for (size_t atom_c = 0; atom_c < numParticle; ++atom_c) {
                        if (atom_r == atom_c)
                            continue;
                        temp += result(r, atom_c * Dim + dim_c);
                    }
                    result[r, c] = -T1(temp);
                }
            }
        }

        if constexpr (ReverseDiff<T>) {
            std::ignore = co_yield std::move(result);
            noImpl(__func__);
        }
        else
            co_return std::move(result);
    }

    template<class Derived>
    auto PairModel<Derived>::virial(const MDCellType& cell) const -> CoDiff<LatticeMatrix> {
        return virial(cell.getLattice(), cell.getPos());
    }
    /**
     * Reference:
     * [1] Chem. Phys. Lett. 421, 138 (2006); https://doi.org/10.1016/J.CPLETT.2006.01.087
     */
    template<class Derived>
    auto PairModel<Derived>::virial(const LatticeMatrix& lattice, const PositionMatrix& pos) const -> CoDiff<LatticeMatrix> {        
        LatticeMatrix result(Dim, Dim, 0);
        auto kernel = [this, &result](size_t i, size_t j, Vec3D r, const T& norm1, const T& norm2) {
            const T f_norm = force_functor(i, j, norm1, norm2);
            const CoDiff<Vec3D> f = r * (f_norm / norm1);
            result += f * r.transpose();
        };
        forPairInCutoff(lattice, pos, kernel);
        result *= reciprocal(MDCellType::getVolume(lattice));
        if constexpr (ReverseDiff<T>) {
            std::ignore = co_yield std::move(result);
            noImpl(__func__);
        }
        else
            co_return std::move(result);
    }

    template<class Derived>
    void PairModel<Derived>::swap(PairModel& __restrict pair) noexcept {
        assert(this != &pair && "[Error]: Self swap is likely a bug");
        cutoff.swap(pair.cutoff);
        squared_cutoff.swap(pair.squared_cutoff);
        pot_shift.swap(pair.pot_shift);
    }

    template<class Derived>
    void PairModel<Derived>::setCutoff(Tv cutoff_) {
        cutoff = std::move(cutoff_);
        squared_cutoff = square(cutoff);
        if constexpr (!IsPotDependOnAtomIndex) {
            constexpr int unused = 0;
            pot_shift = pot_functor(unused, unused, cutoff, squared_cutoff).value();
        }
    }

    template<class Derived>
    void PairModel<Derived>::forPairInCutoff(
            const LatticeMatrix& lattice, const PositionMatrix& cartesianPos, std::invocable<int, int, Vec3D, T, T> auto fn) const {
        const auto& pos = cartesianPos;
        if constexpr (IsSmallCell) {
            const auto range = MDCellType::estimateRange(lattice, cutoff);
            const size_t numParticle = pos.getRow();
            MDCellType::forCellInRange(range, lattice,
                [this, pos, numParticle, &fn](Vec3D delta) {
                    Vec3D r, from;
                    for (size_t i = 0; i < numParticle; ++i) {
                        from = pos.row(i) + delta;
                        for (size_t j = i; j < numParticle; ++j) {
                            auto to = pos.row(j);
                            r = to - from;
                            const T norm2 = r.squaredNorm();
                            const bool isNotSelf = !norm2.isSubNormal();
                            if (isNotSelf && norm2 < squared_cutoff) {
                                const T norm1 = sqrt(norm2);
                                fn(i, j, r, norm1, norm2);
                            }
                        }
                    }
                });
        }
        else {
            const CellListType cellList(lattice, pos, cutoff);
            Array<size_t> arr1{};
            cellList.forCellInList([this, pos, &arr1, &fn, &cellList](Index3D center) {
                cellList.forAtomInCell(center, [&arr1](size_t atom) {
                    arr1.append(atom);
                });
                /* Current cell */ {
                    const size_t length = arr1.getLength();
                    for (size_t i = 0; i + 1 < length; ++i) {
                        const size_t atom1 = arr1[i];
                        const auto from = pos.row(atom1);
                        for (size_t j = i + 1; j < length; ++j) {
                            const size_t atom2 = arr1[j];
                            const auto to = pos.row(atom2);
                            CoDiff<Vec3D> r = to - from;
                            const T norm2 = r.squaredNorm();
                            if (norm2 < squared_cutoff) {
                                const T norm1 = sqrt(norm2);
                                fn(atom1, atom2, r, norm1, norm2);
                            }
                        }
                    }
                }
                Array<size_t> arr2{};
                cellList.forReducedNeighInRange(center, [this, pos, &arr1, &arr2, &fn, &cellList](Vec3D translate, Index3D neigh) {
                    cellList.forAtomInCell(neigh, [&arr2](size_t atom) {
                        arr2.append(atom);
                    });
                    std::sort(arr2.begin(), arr2.end());

                    for (const size_t atom1 : arr1) {
                        const CoDiff<Vec3D> from = pos.row(atom1) - translate;
                        Vec3D r, f(Dim, 0);
                        for (const size_t atom2 : arr2) {
                            const auto to = pos.row(atom2);
                            r = to - from;
                            const T norm2 = r.squaredNorm();
                            if (norm2 < squared_cutoff) {
                                const T norm1 = sqrt(norm2);
                                fn(atom1, atom2, r, norm1, norm2);
                            }
                        }
                    }
                    arr2.clear();
                });
                arr1.clear();
            });
        }
    }

    template<class Derived>
    auto PairModel<Derived>::forceConstImpl(const MDCellType& cell, size_t dof1, size_t dof2) const -> CoDiff<T> {
        static_assert(!IsPotDependOnAtomIndex, "[Error]: It is assumed force not depends on atom index");
        constexpr int unused = 0;
        const size_t atom1 = dof1 / 3U;
        const size_t atom2 = dof2 / 3U;
        assert(atom1 != atom2 && "[Error]: The function is not responsible for this case");
        const size_t dir1 = dof1 % 3U;
        const size_t dir2 = dof2 % 3U;
        const CoDiff<Vec3D> delta = cell.getPos().row(atom2) - cell.getPos().row(atom1);
        const auto squaredNorm = delta.squaredNorm();
        const auto norm = sqrt(squaredNorm);
        const auto factor = delta[dir1] * delta[dir2] / squaredNorm;
        const auto term1 = forceConst_functor(norm, squaredNorm) * factor;
        const auto term2 = force_functor(unused, unused, norm, squaredNorm) / norm * (T(dir1 == dir2 ? 1.0 : 0.0) - factor);
        const auto result = term1 - term2;
        if constexpr (ReverseDiff<T>) {
            auto& r = co_yield result.value();
            result.reverse(r.grad());
        }
        else
            co_return std::move(result);
    }
}

namespace Physica {
    template<class T>
    class Traits<PairModel<T>> {
    public:
        using Derived = T;

        constexpr static bool IsContractable = false;
    };
}
