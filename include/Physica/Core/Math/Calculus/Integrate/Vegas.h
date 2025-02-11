/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "IntegrateImpl/AdaptiveBase.h"

namespace Physica {
    /**
     * Adaptive monte-carlo for high dimensional integration
     *
     * \tparam TakeLn targets divergent integrals and computes logarithmic observations.
     * 
     * Reference:
     * [1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
     * [2] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:414-416
     */
    template<Scalar T, bool TakeLn>
    class Vegas : public AdaptiveBase<T, TakeLn> {
        using This = Vegas<T, TakeLn>;
        using Base = AdaptiveBase<T, TakeLn>;
    protected:
        using Tv = T::ValueType;
        using Trv = Tv::RealType;
        using LossMatrix = DenseMatrix<Trv>;
        using CountArray = Array<Array<int>>;
    private:
        DenseMatrix<Trv> pointGrid;
        Trv compressRate;
        Trv mixBeta;
    protected:
        using Base::from;
        using Base::to;
        using Base::means;
        using Base::vars;
        using Base::loss;

        LossMatrix lossMat;
        CountArray counts;
    public:
        Vegas() = default;
        Vegas(VectorND<Trv> from,
              VectorND<Trv> to,
              int numRefine,
              int numSample,
              int numPoint = 1000,
              Trv compressRate_ = 1.5,
              Trv mixBeta_ = 1);
        Vegas(const This&) = default;
        Vegas(This&&) noexcept = default;
        ~Vegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void mesh_uniform(size_t dim);
        void mesh_tanh(size_t dim, Trv range);

        template<class Functor, RandomGenerator R, class Executor = SeqExecutor>
        Trv warmup(Functor func, int numWarm);
        template<class Functor, RandomGenerator R, class Executor = SeqExecutor>
        void integral(Functor func);
        template<class Functor, RandomGenerator R, class Executor = SeqExecutor>
        [[nodiscard]] Trv calcGridLoss(Functor func) const;

        const H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getDim;
        [[nodiscard]] const auto& getPointGrid() const noexcept { return pointGrid; }
        [[nodiscard]] size_t getNumPoint() const noexcept { return pointGrid.getRow(); }
        /* Setters */
        void setMixBeta(Trv beta) { mixBeta = beta; }
    private:
        void pre_trial();
        template<class Executor>
        void refineGrid();
        Trv calcGridLossImpl() const;
        Trv compress(VectorND<Trv>& vars);

        template<RandomGenerator R>
        [[nodiscard]] std::pair<VectorND<Trv>, VectorND<Trv>> sample(const int* indexes) const;
        template<class Functor, RandomGenerator R, class Executor>
        void trial_normal(Functor func, T& mean, T& var);
        template<class Functor, RandomGenerator R, class Executor>
        void trial_ln(Functor func, T& mean, T& var);
    };
}

#include "IntegrateImpl/VegasImpl.h"
