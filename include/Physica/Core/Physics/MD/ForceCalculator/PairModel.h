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
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType, class PairFunctor>
    class PairModel final {
        using ResultType = typename std::invoke_result<PairFunctor, ScalarType>::type;
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        static_assert(std::is_same<ScalarType, ResultType>::value, "[Error]: Invalid PairFunctor");
    public:
        constexpr static unsigned int Dim = 3;
    private:
        ScalarType cutoff;
        ScalarType squared_cutoff;
        PairFunctor force_functor;
        PairFunctor pot_functor;
    public:
        PairModel(ScalarType cutoff_, PairFunctor functor_, PairFunctor pot_functor_);
        [[nodiscard]] Vector<ScalarType> operator()(MDCell cell) const;
        [[nodiscard]] ScalarType potentialEnergy(MDCell cell) const;
    private:
        Utils::Array<ssize_t, 3> estimateRange(const MDCell& cell) const;
    };

    template<class ScalarType, class PairFunctor>
    PairModel<ScalarType, PairFunctor>::PairModel(ScalarType cutoff_, PairFunctor functor_, PairFunctor pot_functor_)
            : cutoff(std::move(cutoff_))
            , force_functor(std::move(functor_))
            , pot_functor(std::move(pot_functor_)) {
        squared_cutoff = square(cutoff);
    }

    template<class ScalarType, class PairFunctor>
    Vector<ScalarType> PairModel<ScalarType, PairFunctor>::operator()(MDCell cell) const {
        using VectorType = Vector<ScalarType, Dim>;

        const auto& lattice = cell.getLattice();
        const auto& pos = cell.getPos();
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);
        const auto range = estimateRange(cell);
        const size_t numParticle = cell.getNumParticle();

        Vector<ScalarType> force(Dim * numParticle, 0);
        VectorType v1, v2, v3, r;
        for (size_t i = 0; i < numParticle; ++i) {
            auto force_i = force.segment(3 * i, 3 * i + 3);
            auto center = pos.row(i);
            for (size_t j = i; j < numParticle; ++j) {
                auto force_j = force.segment(3 * j, 3 * j + 3);
                auto v = pos.row(j);
                for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                    v1 = v + ScalarType(x) * a1.asVector();
                    for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                        v2 = v1 + ScalarType(y) * a2.asVector();
                        for (ssize_t z = -range[2]; z <= range[2]; ++z) {
                            v3 = v2 + ScalarType(z) * a3.asVector();
                            r = v3 - center;
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
                }
            }
        }
        return force;
    }

    template<class ScalarType, class PairFunctor>
    ScalarType PairModel<ScalarType, PairFunctor>::potentialEnergy(MDCell cell) const {
        using VectorType = Vector<ScalarType, Dim>;

        const auto& lattice = cell.getLattice();
        const auto& pos = cell.getPos();
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);
        const auto range = estimateRange(cell);
        const size_t numParticle = cell.getNumParticle();

        ScalarType result = 0;
        VectorType v1, v2, v3, r;
        for (size_t i = 0; i < numParticle; ++i) {
            auto center = pos.row(i);
            for (size_t j = i; j < numParticle; ++j) {
                auto v = pos.row(j);
                for (ssize_t x = -range[0]; x <= range[0]; ++x) {
                    v1 = v + ScalarType(x) * a1.asVector();
                    for (ssize_t y = -range[1]; y <= range[1]; ++y) {
                        v2 = v1 + ScalarType(y) * a2.asVector();
                        for (ssize_t z = -range[2]; z <= range[2]; ++z) {
                            v3 = v2 + ScalarType(z) * a3.asVector();
                            r = v3 - center;
                            const ScalarType r2 = r.squaredNorm();
                            const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                            if (isNotSelf && r2 < squared_cutoff) {
                                const ScalarType dist = sqrt(r2);
                                result += pot_functor(dist);
                            }
                        }
                    }
                }
            }
        }
        return result;
    }

    template<class ScalarType, class PairFunctor>
    Utils::Array<ssize_t, 3> PairModel<ScalarType, PairFunctor>::estimateRange(const MDCell& cell) const {
        ssize_t max_x, max_y, max_z;
        /* Estimate range */ {
            const ReciprocalCell reciprocal = cell.reciprocal();
            const auto& lattice = reciprocal.getLattice();
            const ScalarType factor = cutoff * (1 / (2 * M_PI));
            max_x = static_cast<ssize_t>(double(factor * lattice.row(0).norm()));
            max_y = static_cast<ssize_t>(double(factor * lattice.row(1).norm()));
            max_z = static_cast<ssize_t>(double(factor * lattice.row(2).norm()));
        }
        return Utils::Array<ssize_t, 3>{max_x, max_y, max_z};
    }
}
