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
        constexpr static unsigned int Dim = MDCellType::Dim;
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
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
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
    Vector<ScalarType> PairModel<ScalarType, PosScalarType, PairFunctor>::force(const MDCellType& cell) const {
        using VectorType = Vector<PosScalarType, Dim>;
        const auto& lattice = cell.getLattice();
        const auto& pos = cell.getPos();
        const auto range = MDCellType::estimateRange(lattice, cutoff);
        const size_t numParticle = cell.getNumParticle();

        Vector<ScalarType> force(Dim * numParticle, 0);
        MDCellType::forCellInRange(range, lattice,
            [this, pos, numParticle, &force](VectorType delta) {
                VectorType r, from;
                for (size_t i = 0; i < numParticle; ++i) {
                    auto force_i = force.segment(3 * i, 3 * i + 3);
                    from = pos.row(i) + delta;
                    for (size_t j = i; j < numParticle; ++j) {
                        auto force_j = force.segment(3 * j, 3 * j + 3);
                        auto to = pos.row(j);
                        r = to.asVector() - from;
                        const ScalarType r2 = r.squaredNorm();
                        const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                        if (isNotSelf && r2 < squared_cutoff) {
                            const ScalarType dist = sqrt(r2);
                            const ScalarType f_norm = force_functor(dist);
                            r *= f_norm / dist;
                            const VectorType& f = r;
                            force_i -= f;
                            force_j += f;
                        }
                    }
                }
            });
        return force;
    }

    template<class ScalarType, class PosScalarType, class PairFunctor>
    ScalarType PairModel<ScalarType, PosScalarType, PairFunctor>::potentialEnergy(const MDCellType& cell) const {
        using VectorType = Vector<PosScalarType, Dim>;
        const auto& pos = cell.getPos();
        const auto range = MDCellType::estimateRange(cell.getLattice(), cutoff);
        const size_t numParticle = cell.getNumParticle();

        ScalarType result = 0;
        MDCellType::forCellInRange(range, cell.getLattice(),
            [this, numParticle, &pos, &result](VectorType delta) {
                ScalarType temp = 0;
                for (size_t i = 0; i < numParticle; ++i) {
                    VectorType from = pos.row(i) + delta;
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
    void PairModel<ScalarType, PosScalarType, PairFunctor>::swap(PairModel& pair) noexcept {
        cutoff.swap(pair.cutoff);
        squared_cutoff.swap(squared_cutoff);
        pot_shift.swap(pot_shift);
        std::swap(force_functor, force_functor);
        std::swap(pot_functor, pot_functor);
    }
}
