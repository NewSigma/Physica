/*
 * Copyright 2024-2026 Weibo He.
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
     *
     * Publication:
     * [3] APL Comput. Phys. 2, 026108 (2026); https://doi.org/10.1063/5.0320242
     */
    template<Scalar T, bool TakeLn = false>
    class Vegas : public AdaptiveBase<T, TakeLn> {
        using This = Vegas<T, TakeLn>;
        using Base = AdaptiveBase<T, TakeLn>;
    protected:
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tv::RealType;
    public:
        using Cube = DenseMatrix<Trv, MatrixMajor::Row, 2, Dynamic>;
    private:
        Cube cube;
        DenseMatrix<Trv> pointGrid;
        Trv compressRate;
        Trv lr;
        Trv momentum;
    protected:
        using Base::means;
        using Base::vars;
        using Base::loss;

        DenseMatrix<Trv> losses;
        Array2D<int> counts;
        Array2D<int> indices;
    public:
        Vegas() = default;
        Vegas(Cube cube,
              int numRefine,
              int numSample,
              int numPoint = 1000,
              Trv compressRate = 1.5,
              Trv lr = 1,
              Trv momentum = 0);
        Vegas(const This&) = default;
        Vegas(This&&) noexcept = default;
        ~Vegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void mesh_uniform(size_t dim);
        void mesh_tanh(size_t dim, Trv range);

        template<RNG R, ExecutePolicy P = Sequential>
        Trv warmup(std::invocable<VectorND<Trv>> auto fn, int numWarm);
        template<RNG R, ExecutePolicy P = Sequential>
        void integral(std::invocable<VectorND<Trv>> auto fn);
        template<RNG R, ExecutePolicy P = Sequential>
        [[nodiscard]] Trv calcGridLoss(std::invocable<VectorND<Trv>> auto fn) const;

        const H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getPointGrid() const noexcept { return pointGrid; }
        [[nodiscard]] size_t getDim() const noexcept { return cube.getCol(); }
        [[nodiscard]] size_t getNumPoint() const noexcept { return pointGrid.getRow(); }
        /* Setters */
        void setLearnRate(Trv lr_) { lr = lr_; }
    private:
        template<RNG R>
        void init();
        template<ExecutePolicy P>
        void compress();
        template<ExecutePolicy P>
        void refineGrid();
        Trv calcGridLossImpl() const;

        template<RNG R>
        [[nodiscard]] Array<VectorND<Trv>, 2> transform(std::span<int> indices) const;
        [[nodiscard]] Vector2D<T> statistic(const VectorND<T>& samples);
        template<RNG R, ExecutePolicy P>
        Vector2D<T> sample(std::invocable<VectorND<Trv>> auto fn);
        template<RNG R, ExecutePolicy P>
        Vector2D<T> sample_ln(std::invocable<VectorND<Trv>> auto fn);
        /* Static members */
        template<class Sample>
        consteval static void checkIntegrand(std::invocable<Sample> auto& fn) noexcept;
    };
}

#include "IntegrateImpl/VegasImpl.h"
