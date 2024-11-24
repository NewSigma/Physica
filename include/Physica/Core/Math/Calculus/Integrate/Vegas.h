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

        DenseMatrix<T> pointGrid;
        VectorND<T> from;
        VectorND<T> to;
        int numRefine;
        int numSample;
        T compressRate;
    protected:
        VectorND<T> means;
        VectorND<T> vars;
    public:
        Vegas() = default;
        Vegas(VectorND<T> from_, VectorND<T> to_, int numRefine_, int numSample_, int numPoint = 1000, T compressRate_ = 1.5);
        Vegas(const This&) = default;
        Vegas(This&&) noexcept = default;
        ~Vegas() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor, RandomGenerator R, class Executor = SequentialExecutor>
        void integral(Functor func);
        template<class Functor, RandomGenerator R>
        [[nodiscard]] T accessMerit(Functor func);

        [[nodiscard]] T calcMean() const;
        [[nodiscard]] T calcDevia() const;
        [[nodiscard]] T calcVar() const;
        [[nodiscard]] T calcSquaredChi() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getPointGrid() const noexcept { return pointGrid; }
        [[nodiscard]] size_t getNumPoint() const noexcept { return pointGrid.getRow(); }
        [[nodiscard]] size_t getDim() const noexcept { return pointGrid.getCol(); }
        [[nodiscard]] int getNumRefine() const noexcept { return numRefine; }
        [[nodiscard]] int getNumSample() const noexcept { return numSample; }
        /* Setters */
        void setNumRefine(int numRefine_);
    protected:
        template<class Executor>
        void refineGrid(DenseMatrix<T>& varsDevia);
        static T accessMeritImpl(const DenseMatrix<T>& varsDevia);
    private:
        template<class Functor, RandomGenerator R, class Executor>
        void trialIntegral(DenseMatrix<T>& varsDevia, int refine, Functor func);
        T compress(VectorND<T>& vars);
    };

    template<Scalar T>
    Vegas<T>::Vegas(VectorND<T> from_, VectorND<T> to_, int numRefine_, int numSample_, int numPoint, T compressRate_)
            : from(std::move(from_))
            , to(std::move(to_))
            , numSample(numSample_)
            , compressRate(compressRate_) {
        assert(from.getLength() == to.getLength() && "[Error]: Inconsistent dim");
        assert(numSample > 0);
        assert(numPoint > 2 && "[Error]: Invalid point number");
        assert(compressRate.isPositive() && "[Error]: Rate = 0 implies no grid refinement");
        assert(compressRate < T(2) && "[Error]: Rate should be ~ 1");
        pointGrid.resize(numPoint, from.getLength());
        for (size_t i = 0; i < getDim(); ++i)
            pointGrid.col(i) = VectorND<T>::linspace(from[i], to[i], numPoint);
        setNumRefine(numRefine_);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::integral(Functor func) {
        using CallResult = std::invoke_result<Functor, VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        DenseMatrix<T> varsDevia(getNumPoint() - 1, getDim(), 0);
        for (int refine = 0; refine < numRefine; ++refine) {
            trialIntegral<Functor, R, Executor>(varsDevia, refine, func);
            refineGrid<Executor>(varsDevia);
        }
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R>
    T Vegas<T>::accessMerit(Functor func) {
        DenseMatrix<T> varsDevia(getNumPoint() - 1, getDim(), 0);
        trialIntegral<Functor, R, SequentialExecutor>(varsDevia, 0, func);
        return accessMeritImpl(varsDevia);
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
            return 0;
        const T mean1 = calcMean();
        return divide(square(means - mean1), vars).sum() / T(numRefine - 1); // Normalize, refer to [2]
    }

    template<Scalar T>
    void Vegas<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pointGrid.swap(obj.pointGrid);
        from.swap(obj.from);
        to.swap(obj.to);
        std::swap(numRefine, obj.numRefine);
        std::swap(numSample, obj.numSample);
        compressRate.swap(obj.compressRate);

        means.swap(obj.means);
        vars.swap(obj.vars);
    }

    template<Scalar T>
    template<class Functor, RandomGenerator R, class Executor>
    void Vegas<T>::trialIntegral(DenseMatrix<T>& varsDevia, int refine, Functor func) {
        const auto indexes = R::getInstance().random_int(getDim() * numSample, 0, varsDevia.getRow() - 1);
        VectorND<T> samples(numSample);
        Executor::parallel_for([&, this](size_t n) {
            VectorND<T> fromX(getDim());
            VectorND<T> deltaX(getDim());
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                fromX[i] = pointGrid(index, i);
                deltaX[i] = pointGrid(index + 1, i);
            }
            deltaX -= fromX;

            const VectorND<T> x = fromX + hadamard(deltaX, VectorND<T>::template random_uniform<R>(getDim()));
            const T y = func(x);
            const T xy = y * (deltaX * T(varsDevia.getRow())).prod();
            samples[n] = xy;
        }, numSample, Executor::getNumThread()).wait();

        Array<Array<int>> counts(varsDevia.getRow(), getDim(), 0);
        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            toNextVariance(vars[refine], means[refine], n, xy);
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indexes[n * getDim() + i];
                toNextMean(varsDevia(index, i), counts[index][i], square(xy));
                counts[index][i] += 1;
            }
        }
    }

    template<Scalar T>
    void Vegas<T>::setNumRefine(int numRefine_) {
        assert(numRefine > 0);
        numRefine = numRefine_;
        means.resize(numRefine, 0);
        vars.resize(numRefine, 0);
    }

    template<Scalar T>
    template<class Executor>
    void Vegas<T>::refineGrid(DenseMatrix<T>& varsDevia) {
        Executor::parallel_for([this, &varsDevia](size_t dim) {
            const T meanVar = compress(varsDevia.asArray()[dim]);
            const bool noData = meanVar.isZero();
            if (noData) [[unlikely]] // No data in the dimension, usually we should have enough samples to avoid it
                return;

            auto oldPoints = pointGrid.col(dim);
            VectorND<T> newPoints(getNumPoint());
            newPoints[0] = oldPoints[0];
            T temp = 0;
            size_t i = 1;
            for (size_t j = 0; i < newPoints.getLength() - 1; ++i) {
                while (temp < meanVar) {
                    assert(j < varsDevia.getRow() && "[Error]: Unexpected not enough vars, this is likely a bug");
                    temp += varsDevia(j, dim);
                    j += 1;
                }
                temp -= meanVar;
                newPoints[i] = oldPoints[j] - temp * (oldPoints[j] - oldPoints[j - 1]) / varsDevia(j - 1, dim);
            }
            newPoints[i] = oldPoints[i];
            oldPoints = newPoints;
        }, getDim(), Executor::getNumThread()).wait();
    }

    template<Scalar T>
    T Vegas<T>::accessMeritImpl(const DenseMatrix<T>& varsDevia) {
        T maxDevia = 0;
        for (size_t i = 0; i < varsDevia.getCol(); ++i) {
            const auto col = varsDevia.col(i);
            const T prior = mean(col);
            maxDevia = std::max(maxDevia, variance(col, prior) / square(prior));
        }
        return maxDevia;
    }

    template<Scalar T>
    T Vegas<T>::compress(VectorND<T>& vars) {
        const T norm1 = vars.norm1();
        const bool noData = norm1.isZero();
        if (noData)
            return 0;

        const VectorND<T> buffer = vars * reciprocal(norm1); // Normalized values fall into a range that is feasible for compression function.
        const Vector3D<T> kernel{1.0 / 8, 6.0 / 8, 1.0 / 8};
        size_t i = 0;
        vars[0] = T(7.0 / 8) * buffer[0] + T(1.0 / 8) * buffer[1];
        for (; i < vars.getLength() - 2; ++i)
            vars[i + 1] = buffer.template segment<3>(i, i + 3) * kernel;
        vars[i + 1] = T(7.0 / 8) * buffer[i] + T(1.0 / 8) * buffer[i + 1];

        vars = pow(divide(vars - T(1), ln(vars)), compressRate);
        return mean(vars);
    }
}
