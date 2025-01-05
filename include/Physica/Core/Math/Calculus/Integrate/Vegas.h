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
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::Core {
    /**
     * Adaptive monte-carlo for high dimensional integration
     * 
     * Reference:
     * [1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
     * [2] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:414-416
     */
    template<Scalar T>
    class Vegas {
        using This = Vegas<T>;
    protected:
        using Tv = T::ValueType;
        using RealValue = Tv::RealType;
        using LossMatrix = DenseMatrix<RealValue>;
        using CountArray = Array<Array<int>>;
    private:
        DenseMatrix<RealValue> pointGrid;
        VectorND<RealValue> from;
        VectorND<RealValue> to;
        int numRefine;
        int numSample;
        RealValue compressRate;
        RealValue mixBeta;
    protected:
        LossMatrix lossMat;
        CountArray counts;

        VectorND<T> means;
        VectorND<T> vars;
        VectorND<RealValue> loss;
    public:
        Vegas() = default;
        Vegas(VectorND<RealValue> from_,
              VectorND<RealValue> to_,
              int numRefine_,
              int numSample_,
              int numPoint = 1000,
              RealValue compressRate_ = 1.5,
              RealValue mixBeta_ = 1);
        Vegas(const This&) = default;
        Vegas(This&&) noexcept = default;
        ~Vegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void mesh_uniform(size_t dim);
        void mesh_tanh(size_t dim, RealValue range);

        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void warmup(Functor func, int numWarm);
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void integral(Functor func);
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        [[nodiscard]] RealValue calcGridLoss(Functor func) const;

        [[nodiscard]] T calcMean() const;
        [[nodiscard]] T calcDevia() const;
        [[nodiscard]] T calcVar() const;
        [[nodiscard]] T calcSquaredChi() const;
    #ifdef PHYSICA_HDF5
        const H5Group read(const H5Location& loc, const char* name);
        H5Group write(H5Location& loc, const char* name) const;
    #endif
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getPointGrid() const noexcept { return pointGrid; }
        [[nodiscard]] size_t getNumPoint() const noexcept { return pointGrid.getRow(); }
        [[nodiscard]] size_t getDim() const noexcept { return pointGrid.getCol(); }
        [[nodiscard]] int getNumRefine() const noexcept { return numRefine; }
        [[nodiscard]] int getNumSample() const noexcept { return numSample; }
        [[nodiscard]] const VectorND<T>& getMeans() const noexcept { return means; }
        [[nodiscard]] const VectorND<T>& getVars() const noexcept { return vars; }
        [[nodiscard]] const VectorND<RealValue>& getLoss() const noexcept { return loss; }
        /* Setters */
        void setNumRefine(int numRefine_);
        void setMixBeta(RealValue beta) { mixBeta = beta; }
    protected:
        void pre_trial();
        template<RandomGenerator R>
        [[nodiscard]] std::pair<VectorND<RealValue>, VectorND<RealValue>> sample(const int* indexes) const;
        template<class Executor>
        void refineGrid();

        RealValue calcGridLossImpl() const;
    private:
        template<class Functor, RandomGenerator R, class Executor>
        std::pair<T, T> trialIntegral(Functor func);
        RealValue compress(VectorND<RealValue>& vars);
    };

    template<Scalar T>
    Vegas<T>::Vegas(VectorND<RealValue> from_, VectorND<RealValue> to_, int numRefine_, int numSample_, int numPoint, RealValue compressRate_, RealValue mixBeta_)
            : from(std::move(from_))
            , to(std::move(to_))
            , numSample(numSample_)
            , compressRate(compressRate_)
            , mixBeta(mixBeta_) {
        assert(from.getLength() == to.getLength() && "[Error]: Inconsistent dim");
        assert(numSample > 0);
        assert(numPoint > 2 && "[Error]: Invalid point number");
        assert(compressRate.isPositive() && "[Error]: Rate = 0 implies no grid refinement");
        assert(compressRate < RealValue(2) && "[Error]: Rate should be ~ 1");
        assert(mixBeta.isPositive() && (mixBeta <= RealValue(1)) && "[Error]: Invalid beta");
        pointGrid.resize(numPoint, from.getLength());
        lossMat.resize(numPoint - 1, getDim());
        counts.resize(numPoint - 1, getDim());
        setNumRefine(numRefine_);
        for (size_t i = 0; i < getDim(); ++i)
            mesh_uniform(i);
    }

    template<Scalar T>
    void Vegas<T>::mesh_uniform(size_t dim) {
        pointGrid.col(dim) = VectorND<RealValue>::linspace(from[dim], to[dim], getNumPoint());
    }

    template<Scalar T>
    void Vegas<T>::mesh_tanh(size_t dim, RealValue range) {
        assert(range.isPositive());
        const auto k = (to[dim] - from[dim]) * 0.5;
        const auto b = (to[dim] + from[dim]) * 0.5;
        auto p = VectorND<RealValue>::linspace(-range, range, getNumPoint());
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

        for (int _ = 0; _ < numWarm; ++_) {
            pre_trial();
            trialIntegral<Functor, R, Executor>(func);
            refineGrid<Executor>();
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::integral(Functor func) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        for (int refine = 0; refine < numRefine; ++refine) {
            pre_trial();
            auto pair = trialIntegral<Functor, R, Executor>(func);
            means[refine] = std::move(pair.first);
            vars[refine] = std::move(pair.second);
            loss[refine] = calcGridLossImpl();
            refineGrid<Executor>();
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    Vegas<T>::RealValue Vegas<T>::calcGridLoss(Functor func) const {
        pre_trial();
        trialIntegral<Functor, R, Executor>(func);
        return calcGridLossImpl();
    }

    template<Scalar T>
    T Vegas<T>::calcMean() const {
        return reciprocal(vars) * means / reciprocal(vars).sum();
    }

    template<Scalar T>
    T Vegas<T>::calcVar() const {
        return reciprocal(reciprocal(vars).sum());
    }

    template<Scalar T>
    T Vegas<T>::calcDevia() const {
        return sqrt(calcVar());
    }

    template<Scalar T>
    T Vegas<T>::calcSquaredChi() const {
        if (numRefine == 1)
            return 1;
        const T mean1 = calcMean();
        return (RealValue(numSample) / RealValue(numRefine - 1)) * divide((means - mean1).squaredNorms(), vars).sum(); // Normalize, refer to [2]
    }
#ifdef PHYSICA_HDF5
    template<Scalar T>
    const H5Group Vegas<T>::read(const H5Location& loc, const char* name) {
        const auto group = loc.openGroup(name);
        pointGrid.read(group, "Grid");
        from.read(group, "From");
        to.read(group, "To");

        group.readAttr("NumRefine", numRefine);
        group.readAttr("NumSample", numSample);
        group.readAttr("CompressRate", compressRate);
        group.readAttr("MixBeta", mixBeta);

        lossMat.resize(getNumPoint() - 1, getDim());
        counts.resize(getNumPoint() - 1, getDim());

        means.read(group, "Means");
        vars.read(group, "Vars");
        loss.read(group, "Loss");
        return group;
    }

    template<Scalar T>
    H5Group Vegas<T>::write(H5Location& loc, const char* name) const {
        auto group = loc.openGroup(name);
        pointGrid.write(group, "Grid");
        from.write(group, "From");
        to.write(group, "To");

        group.writeAttr("NumRefine", numRefine);
        group.writeAttr("NumSample", numSample);
        group.writeAttr("CompressRate", compressRate);
        group.writeAttr("MixBeta", mixBeta);

        means.write(group, "Means");
        vars.write(group, "Vars");
        loss.write(group, "Loss");
        return group;
    }
#endif
    template<Scalar T>
    void Vegas<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pointGrid.swap(obj.pointGrid);
        from.swap(obj.from);
        to.swap(obj.to);
        std::swap(numRefine, obj.numRefine);
        std::swap(numSample, obj.numSample);
        compressRate.swap(obj.compressRate);
        mixBeta.swap(obj.mixBeta);

        lossMat.swap(obj.lossMat);
        counts.swap(obj.counts);

        means.swap(obj.means);
        vars.swap(obj.vars);
        loss.swap(obj.loss);
    }

    template<Scalar T>
    void Vegas<T>::setNumRefine(int numRefine_) {
        assert(numRefine_ > 0);
        numRefine = numRefine_;
        means.resize(numRefine);
        vars.resize(numRefine);
        loss.resize(numRefine);
        means = RealValue(0);
        vars = RealValue(0);
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
    auto Vegas<T>::sample(const int* indexes) const -> std::pair<VectorND<RealValue>, VectorND<RealValue>> {
        VectorND<RealValue> x(getDim());
        VectorND<RealValue> deltas(getDim());
        for (size_t i = 0; i < getDim(); ++i) {
            int index = indexes[i];
            x[i] = pointGrid(index, i);
            deltas[i] = pointGrid(index + 1, i);
        }
        deltas -= x;
        x += hadamard(deltas, VectorND<RealValue>::template random_uniform<R>(getDim()));
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
            VectorND<RealValue> newPoints(getNumPoint());
            newPoints[0] = oldPoints[0];
            RealValue temp = 0;
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
            oldPoints = newPoints * mixBeta + oldPoints * (RealValue(1) - mixBeta);
        }, getDim(), Executor::getNumThread()).wait();
    }

    template<Scalar T>
    auto Vegas<T>::calcGridLossImpl() const -> RealValue {
        RealValue maxVar = 0;
        for (size_t i = 0; i < lossMat.getCol(); ++i) {
            const auto col = lossMat.col(i);
            const RealValue prior = mean(col);
            maxVar = std::max(maxVar, variance(col, prior) / prior.squaredNorm());
        }
        return sqrt(maxVar);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    std::pair<T, T> Vegas<T>::trialIntegral(Functor func) {
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, getNumPoint() - 2);
        VectorND<T> samples(numSample);
        Executor::parallel_for([&, this](size_t n) {
            const auto pair = sample<R>(indexes.data_ptr(n * getDim()));
            const auto& x = pair.first;
            const auto& deltas = pair.second;
            const T y = func(x);
            const T xy = y * (deltas * RealValue(getNumPoint())).prod();
            samples[n] = xy;
        }, numSample, Executor::getNumThread()).wait();

        T mean = 0, var = 0;
        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            toNextVariance(var, mean, n, xy);
            // Loss has minimal value to avoid the grid size reducing to 0
            const auto l = std::max(xy.value().squaredNorm(), RealValue(std::numeric_limits<T>::min()));
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                toNextMean(lossMat(index, i), counts[index][i], l);
                counts[index][i] += 1;
            }
        }
        return std::make_pair(std::move(mean), std::move(var));
    }

    template<Scalar T>
    Vegas<T>::RealValue Vegas<T>::compress(VectorND<RealValue>& vars) {
        const auto sum = vars.sum();
        const bool noData = sum.isZero();
        if (noData)
            return 0;

        const VectorND<RealValue> buffer = vars * reciprocal(sum); // Normalized values fall into a range that is feasible for compression function.
        const Vector3D<RealValue> kernel{1.0 / 8, 6.0 / 8, 1.0 / 8};
        size_t i = 0;
        vars[0] = RealValue(7.0 / 8) * buffer[0] + RealValue(1.0 / 8) * buffer[1];
        for (; i < vars.getLength() - 2; ++i)
            vars[i + 1] = buffer.template segment<3>(i, i + 3) * kernel;
        vars[i + 1] = RealValue(1.0 / 8) * buffer[i] + RealValue(7.0 / 8) * buffer[i + 1];

        vars = pow(divide(vars - RealValue(1), ln(vars)), compressRate);
        return mean(vars);
    }
}
