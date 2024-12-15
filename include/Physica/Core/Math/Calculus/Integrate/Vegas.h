/*
 * Copyright 2024 Weibo He.
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
        using ValueType = T::ValueType;
        using RealValue = ValueType::RealType;
        using LossMatrix = DenseMatrix<RealValue>;
    private:
        DenseMatrix<RealValue> pointGrid;
        VectorND<RealValue> from;
        VectorND<RealValue> to;
        int numRefine;
        int numSample;
        RealValue compressRate;
        RealValue mixBeta;
    protected:
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
              RealValue mixBeta_ = 0);
        Vegas(const This&) = default;
        Vegas(This&&) noexcept = default;
        ~Vegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
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
    protected:
        template<class Executor>
        void refineGrid(LossMatrix& lossMat);
        static RealValue calcGridLossImpl(const LossMatrix& lossMat);
    private:
        template<class Functor, RandomGenerator R, class Executor>
        std::pair<T, T> trialIntegral(LossMatrix& lossMat, Functor func) const;
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
        assert(!mixBeta.isNegative() && (mixBeta < RealValue(1)) && "[Error]: Invalid beta");
        pointGrid.resize(numPoint, from.getLength());
        for (size_t i = 0; i < getDim(); ++i)
            pointGrid.col(i) = VectorND<RealValue>::linspace(from[i], to[i], numPoint);
        setNumRefine(numRefine_);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::warmup(Functor func, int numWarm) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");
        assert(numWarm >= 0 && "[Error]: Invalid param");

        LossMatrix lossMat(getNumPoint() - 1, getDim());
        for (int _ = 0; _ < numWarm; ++_) {
            auto pair = trialIntegral<Functor, R, Executor>(lossMat, func);
            refineGrid<Executor>(lossMat);
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::integral(Functor func) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        LossMatrix lossMat(getNumPoint() - 1, getDim());
        for (int refine = 0; refine < numRefine; ++refine) {
            auto pair = trialIntegral<Functor, R, Executor>(lossMat, func);
            means[refine] = std::move(pair.first);
            vars[refine] = std::move(pair.second);
            loss[refine] = calcGridLossImpl(lossMat);
            refineGrid<Executor>(lossMat);
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    Vegas<T>::RealValue Vegas<T>::calcGridLoss(Functor func) const {
        LossMatrix lossMat(getNumPoint() - 1, getDim());
        trialIntegral<Functor, R, Executor>(lossMat, func);
        return calcGridLossImpl(lossMat);
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
        return (RealValue(numSample) / RealValue(numRefine - 1)) * divide(toSquaredNormVector(means - mean1), vars).sum(); // Normalize, refer to [2]
    }
#ifdef PHYSICA_HDF5
    template<Scalar T>
    const H5Group Vegas<T>::read(const H5Location& loc, const char* name) {
        const auto group = loc.openGroup(name);
        pointGrid.read(group, "Grid");
        from.read(group, "From");
        to.read(group, "To");

        means.read(group, "Means");
        vars.read(group, "Vars");
        loss.read(group, "Loss");

        group.readAttr("NumRefine", numRefine);
        group.readAttr("NumSample", numSample);
        group.readAttr("CompressRate", compressRate);
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

        means.swap(obj.means);
        vars.swap(obj.vars);
        loss.swap(obj.loss);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    std::pair<T, T> Vegas<T>::trialIntegral(LossMatrix& lossMat, Functor func) const {
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, lossMat.getRow() - 1);
        VectorND<T> samples(numSample);
        Executor::parallel_for([&, this](size_t n) {
            VectorND<RealValue> fromX(getDim());
            VectorND<RealValue> deltaX(getDim());
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                fromX[i] = pointGrid(index, i);
                deltaX[i] = pointGrid(index + 1, i);
            }
            deltaX -= fromX;

            const VectorND<RealValue> x = fromX + hadamard(deltaX, VectorND<RealValue>::template random_uniform<R>(getDim()));
            const T y = func(x);
            const T xy = y * (deltaX * RealValue(lossMat.getRow())).prod();
            samples[n] = xy;
        }, numSample, Executor::getNumThread()).wait();

        Array<Array<int>> counts(lossMat.getRow(), getDim(), 0);
        lossMat = std::numeric_limits<T>::min(); // Initial value avoids situation where number of samples is too small to sample effective data

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
    template<class Executor>
    void Vegas<T>::refineGrid(LossMatrix& lossMat) {
        Executor::parallel_for([this, &lossMat](size_t dim) {
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
                    assert(j < lossMat.getRow() && "[Error]: Unexpected not enough vars, this is likely a bug");
                    temp += lossMat(j, dim);
                    j += 1;
                }
                temp -= meanVar;
                newPoints[i] = oldPoints[j] - temp * (oldPoints[j] - oldPoints[j - 1]) / lossMat(j - 1, dim);
            }
            newPoints[i] = oldPoints[i];
            oldPoints = newPoints * (RealValue(1) - mixBeta) + oldPoints * mixBeta;
        }, getDim(), Executor::getNumThread()).wait();
    }

    template<Scalar T>
    Vegas<T>::RealValue Vegas<T>::calcGridLossImpl(const LossMatrix& lossMat) {
        RealValue maxVar = 0;
        for (size_t i = 0; i < lossMat.getCol(); ++i) {
            const auto col = lossMat.col(i);
            const RealValue prior = mean(col);
            maxVar = std::max(maxVar, variance(col, prior) / prior.squaredNorm());
        }
        return sqrt(maxVar);
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
