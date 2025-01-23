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

namespace Physica::Core {
    /**
     * Adaptive monte-carlo for high dimensional integration
     * 
     * Reference:
     * [1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
     * [2] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:414-416
     */
    template<Scalar T>
    class Vegas : public AdaptiveBase<T> {
        using This = Vegas<T>;
        using Base = AdaptiveBase<T>;
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

        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void warmup(Functor func, int numWarm);
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void integral(Functor func);
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        [[nodiscard]] Trv calcGridLoss(Functor func) const;
    #ifdef PHYSICA_HDF5
        const H5Group read(const H5Location& loc, const char* name);
        H5Group write(H5Location& loc, const char* name) const;
    #endif
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getDim;
        [[nodiscard]] const auto& getPointGrid() const noexcept { return pointGrid; }
        [[nodiscard]] size_t getNumPoint() const noexcept { return pointGrid.getRow(); }
        /* Setters */
        void setMixBeta(Trv beta) { mixBeta = beta; }
    protected:
        void pre_trial();
        template<RandomGenerator R>
        [[nodiscard]] std::pair<VectorND<Trv>, VectorND<Trv>> sample(const int* indexes) const;
        template<class Executor>
        void refineGrid();

        Trv calcGridLossImpl() const;
    private:
        template<class Functor, RandomGenerator R, class Executor>
        void trialIntegral(Functor func, T& mean, T& var);
        Trv compress(VectorND<Trv>& vars);
    };

    template<Scalar T>
    Vegas<T>::Vegas(VectorND<Trv> from, VectorND<Trv> to, int numRefine, int numSample, int numPoint, Trv compressRate_, Trv mixBeta_)
            : Base(std::move(from), std::move(to), numRefine, numSample)
            , compressRate(compressRate_)
            , mixBeta(mixBeta_) {
        assert(numPoint > 2 && "[Error]: Invalid point number");
        assert(compressRate.isPositive() && "[Error]: Rate = 0 implies no grid refinement");
        assert(compressRate < Trv(2) && "[Error]: Rate should be ~ 1");
        assert(mixBeta.isPositive() && (mixBeta <= Trv(1)) && "[Error]: Invalid beta");
        pointGrid.resize(numPoint, getDim());
        lossMat.resize(numPoint - 1, getDim());
        counts.resize(numPoint - 1, getDim());
        
        for (size_t i = 0; i < getDim(); ++i)
            mesh_uniform(i);
    }

    template<Scalar T>
    void Vegas<T>::mesh_uniform(size_t dim) {
        pointGrid.col(dim) = VectorND<Trv>::linspace(from[dim], to[dim], getNumPoint());
    }

    template<Scalar T>
    void Vegas<T>::mesh_tanh(size_t dim, Trv range) {
        assert(range.isPositive());
        const auto k = (to[dim] - from[dim]) * 0.5;
        const auto b = (to[dim] + from[dim]) * 0.5;
        auto p = VectorND<Trv>::linspace(-range, range, getNumPoint());
        p = tanh(p);
        p[0] = -1;
        p[getNumPoint() - 1] = 1;
        pointGrid.col(dim) = k * p + b;
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::warmup(Functor func, int numWarm) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");
        assert(numWarm >= 0 && "[Error]: Invalid param");

        T mean, var;
        for (int _ = 0; _ < numWarm; ++_) {
            pre_trial();
            trialIntegral<Functor, R, Executor>(func, mean, var);
            refineGrid<Executor>();
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::integral(Functor func) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        const int numRefine = Base::getNumRefine();
        for (int refine = 0; refine < numRefine; ++refine) {
            pre_trial();
            trialIntegral<Functor, R, Executor>(func, means[refine], vars[refine]);
            loss[refine] = calcGridLossImpl();
            refineGrid<Executor>();
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    Vegas<T>::Trv Vegas<T>::calcGridLoss(Functor func) const {
        T mean, var;
        pre_trial();
        trialIntegral<Functor, R, Executor>(func, mean, var);
        return calcGridLossImpl();
    }
#ifdef PHYSICA_HDF5
    template<Scalar T>
    const H5Group Vegas<T>::read(const H5Location& loc, const char* name) {
        const auto group = Base::read(loc, name);
        group.readAttr("CompressRate", compressRate);
        group.readAttr("MixBeta", mixBeta);

        pointGrid.read(group, "Grid");
        lossMat.resize(getNumPoint() - 1, getDim());
        counts.resize(getNumPoint() - 1, getDim());
        return group;
    }

    template<Scalar T>
    H5Group Vegas<T>::write(H5Location& loc, const char* name) const {
        auto group = Base::write(loc, name);
        group.writeAttr("CompressRate", compressRate);
        group.writeAttr("MixBeta", mixBeta);

        pointGrid.write(group, "Grid");
        return group;
    }
#endif
    template<Scalar T>
    void Vegas<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        pointGrid.swap(obj.pointGrid);
        compressRate.swap(obj.compressRate);
        mixBeta.swap(obj.mixBeta);
        lossMat.swap(obj.lossMat);
        counts.swap(obj.counts);
    }

    template<Scalar T>
    void Vegas<T>::pre_trial() {
        for (auto& arr : counts)
            for (int& elem : arr)
                elem = 0;
        lossMat = std::numeric_limits<T>::min(); // Initial value avoids situation where number of samples is too small to sample effective data
    }

    template<Scalar T>
    template<RandomGenerator R>
    auto Vegas<T>::sample(const int* indexes) const -> std::pair<VectorND<Trv>, VectorND<Trv>> {
        VectorND<Trv> x(getDim());
        VectorND<Trv> deltas(getDim());
        for (size_t i = 0; i < getDim(); ++i) {
            int index = indexes[i];
            x[i] = pointGrid(index, i);
            deltas[i] = pointGrid(index + 1, i);
        }
        deltas -= x;
        x += hadamard(deltas, VectorND<Trv>::template random_uniform<R>(getDim()));
        return std::make_pair(std::move(x), std::move(deltas));
    }

    template<Scalar T>
    template<class Executor>
    void Vegas<T>::refineGrid() {
        Executor::parallel_for([this](size_t dim) {
            const auto meanVar = compress(lossMat.asArray()[dim]);
            const bool noData = meanVar.isZero();
            if (noData) [[unlikely]] // No data in the dimension, usually we should have enough samples to avoid it
                return;

            auto oldPoints = pointGrid.col(dim);
            VectorND<Trv> newPoints(getNumPoint());
            newPoints[0] = oldPoints[0];
            Trv temp = 0;
            size_t i = 1;
            for (size_t j = 0; i < newPoints.getLength() - 1; ++i) {
                while (temp < meanVar) {
                    assert(j < getNumPoint() && "[Error]: Unexpected not enough vars, this is likely a bug");
                    temp += lossMat(j, dim);
                    j += 1;
                }
                temp -= meanVar;
                newPoints[i] = oldPoints[j] - temp * (oldPoints[j] - oldPoints[j - 1]) / lossMat(j - 1, dim);
            }
            newPoints[i] = oldPoints[i];
            oldPoints = newPoints * mixBeta + oldPoints * (Trv(1) - mixBeta);
        }, getDim(), Executor::getNumThread()).wait();
    }

    template<Scalar T>
    auto Vegas<T>::calcGridLossImpl() const -> Trv {
        Trv sumvar = 0;
        for (size_t i = 0; i < lossMat.getCol(); ++i) {
            const auto col = lossMat.col(i);
            const Trv prior = mean(col);
            sumvar += sqrt(variance(col, prior) / prior.squaredNorm() / Trv(getNumPoint()));
        }
        return sumvar / Trv(getDim());
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::trialIntegral(Functor func, T& mean, T& var) {
        const int numSample = Base::getNumSample();
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, getNumPoint() - 2);
        VectorND<T> samples(numSample);
        Executor::parallel_for([&, this](size_t n) {
            const auto pair = sample<R>(indexes.data_ptr(n * getDim()));
            const auto& x = pair.first;
            const auto& deltas = pair.second;
            const T y = func(x);
            const T xy = y * (deltas * Trv(getNumPoint())).prod();
            samples[n] = xy;
        }, numSample, Executor::getNumThread()).wait();

        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            toNextVariance(var, mean, n, xy);
            // Loss has minimal value to avoid the grid size reducing to 0
            const auto l = std::max(xy.value().squaredNorm(), Trv(std::numeric_limits<T>::min()));
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                toNextMean(lossMat(index, i), counts[index][i], l);
                counts[index][i] += 1;
            }
        }
    }

    template<Scalar T>
    Vegas<T>::Trv Vegas<T>::compress(VectorND<Trv>& vars) {
        const auto sum = vars.sum();
        const bool noData = sum.isZero();
        if (noData)
            return 0;

        const VectorND<Trv> buffer = vars * reciprocal(sum); // Normalized values fall into a range that is feasible for compression function.
        const Vector3D<Trv> kernel{1.0 / 8, 6.0 / 8, 1.0 / 8};
        size_t i = 0;
        vars[0] = Trv(7.0 / 8) * buffer[0] + Trv(1.0 / 8) * buffer[1];
        for (; i < vars.getLength() - 2; ++i)
            vars[i + 1] = buffer.template segment<3>(i, i + 3) * kernel;
        vars[i + 1] = Trv(1.0 / 8) * buffer[i] + Trv(7.0 / 8) * buffer[i + 1];

        vars = pow(divide(vars - Trv(1), ln(vars)), compressRate);
        return mean(vars);
    }
}
