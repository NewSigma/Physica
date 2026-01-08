/*
 * Copyright 2025 Weibo He.
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

#include "../Vegas.h"

namespace Physica {
    template<Scalar T, bool TakeLn>
    Vegas<T, TakeLn>::Vegas(VectorND<Trv> from, VectorND<Trv> to, int numRefine, int numSample, int numPoint, Trv compressRate, Trv mixBeta)
            : Base(std::move(from), std::move(to), numRefine, numSample)
            , pointGrid(numPoint, getDim())
            , compressRate(compressRate)
            , mixBeta(mixBeta)
            , lossMat(numPoint - 1, getDim())
            , counts(numPoint - 1, getDim()) {
        assert(numPoint > 2 && "[Error]: Invalid point number");
        assert(compressRate.isPositive() && "[Error]: Rate = 0 implies no grid refinement");
        assert(compressRate < Trv(2) && "[Error]: Rate should be ~ 1");
        assert(mixBeta.isPositive() && (mixBeta <= Trv(1)) && "[Error]: Invalid beta");
        for (size_t i = 0; i < getDim(); ++i)
            mesh_uniform(i);
    }

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::mesh_uniform(size_t dim) {
        pointGrid.col(dim) = VectorND<Trv>::linspace(from[dim], to[dim], getNumPoint());
    }

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::mesh_tanh(size_t dim, Trv range) {
        assert(range.isPositive());
        const auto k = (to[dim] - from[dim]) * 0.5;
        const auto b = (to[dim] + from[dim]) * 0.5;
        auto p = VectorND<Trv>::linspace(-range, range, getNumPoint());
        p = tanh(p);
        p[0] = -1;
        p[getNumPoint() - 1] = 1;
        pointGrid.col(dim) = k * p + b;
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    auto Vegas<T, TakeLn>::warmup(std::invocable<VectorND<Trv>> auto fn, int numWarm) -> Trv {
        using CallResult = std::invoke_result<decltype(fn), VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");
        assert(numWarm >= 0 && "[Error]: Invalid param");

        T mean, var;
        for (int _ = 0; _ < numWarm; ++_) {
            pre_trial();
            if constexpr (TakeLn)
                trial_ln<R, P>(fn, mean, var);
            else
                trial_normal<R, P>(fn, mean, var);
            refineGrid<P>();
        }
        return calcGridLossImpl();
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    void Vegas<T, TakeLn>::integral(std::invocable<VectorND<Trv>> auto fn) {
        using CallResult = std::invoke_result<decltype(fn), VectorND<T>>::type;
        static_assert(std::is_same<CallResult, T>::value, "[Error]: Invalid functor");

        const int numRefine = Base::getNumRefine();
        T mean = 0, var = 0;
        for (int refine = 0; refine < numRefine; ++refine) {
            pre_trial();
            if constexpr (TakeLn)
                trial_ln<R, P>(fn, mean, var);
            else
                trial_normal<R, P>(fn, mean, var);
            means[refine] = mean;
            vars[refine] = var;
            loss[refine] = calcGridLossImpl();
            refineGrid<P>();
        }
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    auto Vegas<T, TakeLn>::calcGridLoss(std::invocable<VectorND<Trv>> auto fn) const -> Trv {
        T mean, var;
        pre_trial();
        if constexpr (TakeLn)
            trial_ln<R, P>(fn, mean, var);
        else
            trial_normal<R, P>(fn, mean, var);
        return calcGridLossImpl();
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, bool TakeLn>
    const H5Group Vegas<T, TakeLn>::read(const H5Loc& loc, const char* name) {
        const auto group = Base::read(loc, name);
        group.readAttr("CompressRate", compressRate);
        group.readAttr("MixBeta", mixBeta);

        pointGrid.read(group, "Grid");
        lossMat.resize(getNumPoint() - 1, getDim());
        counts.resize(getNumPoint() - 1, getDim());
        return group;
    }

    template<Scalar T, bool TakeLn>
    H5Group Vegas<T, TakeLn>::write(H5Loc& loc, const char* name) const {
        auto group = Base::write(loc, name);
        group.writeAttr("CompressRate", compressRate);
        group.writeAttr("MixBeta", mixBeta);

        pointGrid.write(group, "Grid");
        return group;
    }
#endif

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        pointGrid.swap(obj.pointGrid);
        compressRate.swap(obj.compressRate);
        mixBeta.swap(obj.mixBeta);
        lossMat.swap(obj.lossMat);
        counts.swap(obj.counts);
    }

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::pre_trial() {
        lossMat = Trv(std::numeric_limits<T>::min()); // Initial value avoids situation where number of samples is too small to sample effective data
        counts.zeros();
    }

    template<Scalar T, bool TakeLn>
    auto Vegas<T, TakeLn>::compress(Vector auto&& vars) -> Trv {
        assert(vars.min().isPositive());
        const auto sum = vars.sum();
        const bool noData = sum.isZero();
        if (noData)
            return 0;

        const VectorND<Trv> buffer = vars * reciprocal(sum); // Normalized values fall into a range that is feasible for compression function.
        const Vector3D<Trv> kernel{1.0 / 8, 6.0 / 8, 1.0 / 8};
        size_t i = 0;
        vars[0] = fma(Trv(7.0 / 8), buffer[0], Trv(1.0 / 8) * buffer[1]);
        for (; i < vars.getLength() - 2; ++i)
            vars[i + 1] = buffer.template segment<3>(i, i + 3) * kernel;
        vars[i + 1] = fma(Trv(7.0 / 8), buffer[i + 1], Trv(1.0 / 8) * buffer[i]);

        vars = pow(divide(vars - Trv(1), ln(vars)), compressRate);
        return vars.mean();
    }

    template<Scalar T, bool TakeLn>
    template<ExecutePolicy P>
    void Vegas<T, TakeLn>::refineGrid() {
        parallel_for<P>([this](size_t dim) {
            const auto meanVar = compress(lossMat.col(dim));
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
                    temp += lossMat[j, dim];
                    j += 1;
                }
                temp -= meanVar;
                Trv delta = oldPoints[j] - oldPoints[j - 1];
                newPoints[i] = fma(temp / lossMat[j - 1, dim], -delta, oldPoints[j]);
            }
            newPoints[i] = oldPoints[i];
            oldPoints = newPoints * mixBeta + oldPoints * (Trv(1) - mixBeta);
        }, getDim(), 0).wait();
    }

    template<Scalar T, bool TakeLn>
    auto Vegas<T, TakeLn>::calcGridLossImpl() const -> Trv {
        Trv sumvar = 0;
        for (size_t i = 0; i < lossMat.getCol(); ++i) {
            const auto col = lossMat.col(i);
            const Trv prior = col.mean();
            sumvar += sqrt(col.variance(prior) / prior.squaredNorm() / Trv(getNumPoint()));
        }
        return sumvar / Trv(getDim());
    }

    template<Scalar T, bool TakeLn>
    template<RNG R>
    auto Vegas<T, TakeLn>::sample(const int* indices) const -> std::pair<VectorND<Trv>, VectorND<Trv>> {
        VectorND<Trv> x(getDim());
        VectorND<Trv> deltas(getDim());
        for (size_t i = 0; i < getDim(); ++i) {
            int index = indices[i];
            x[i] = pointGrid[index, i];
            deltas[i] = pointGrid[index + 1, i];
        }
        deltas -= x;
        deltas.clamp_min(std::numeric_limits<T>::min());
        x += hadamard(deltas, VectorND<Trv>::template random_uniform<R>(getDim()));
        return std::make_pair(std::move(x), std::move(deltas));
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    void Vegas<T, TakeLn>::trial_normal(std::invocable<VectorND<Trv>> auto fn, T& mean, T& var) {
        const int numSample = Base::getNumSample();
        const auto indices = R::getInstance().random_int(getDim() * numSample, 0, getNumPoint() - 2);
        VectorND<T> samples(numSample);
        parallel_for<P>([&, this](size_t n) {
            const auto pair = sample<R>(indices.data_ptr(n * getDim()));
            const auto& x = pair.first;
            const auto& deltas = pair.second;
            const T y = fn(x);
            assert(y.isFinite() && "[Error]: Bad value");
            const T xy = y * (deltas * Trv(getNumPoint())).prod();
            samples[n] = xy;
        }, numSample, 0).wait();

        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            var.toNextVariance(mean, n, xy);

            const auto l = xy.value().squaredNorm();
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indices[n * getDim() + i];
                lossMat[index, i].toNextMean(counts[index, i], l);
                counts[index, i] += 1;
            }
        }
        // Loss has minimal value to avoid the grid size reducing to 0
        lossMat.clamp_min(std::numeric_limits<Trv>::min());
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    void Vegas<T, TakeLn>::trial_ln(std::invocable<VectorND<Trv>> auto fn, T& mean, T& var) {
        const int numSample = Base::getNumSample();
        const auto indices = R::getInstance().random_int(getDim() * numSample, 0, getNumPoint() - 2);
        VectorND<T> samples(numSample);
        parallel_for<P>([&, this](size_t n) {
            const auto pair = sample<R>(indices.data_ptr(n * getDim()));
            const auto& x = pair.first;
            const auto& deltas = pair.second;
            const T lny = fn(x);
            assert(lny.isFinite() && "[Error]: Bad value");
            T lnxy = lny + ln(deltas).sum();
            lnxy.value() = fma(Trv(getDim()), ln(Trv(getNumPoint())), lnxy.value());
            samples[n] = lnxy;
        }, numSample, 0).wait();

        Tv maxSample;
        if constexpr (T::isComplex)
            maxSample = samples.reals().max().value();
        else
            maxSample = samples.max().value(); // Assuming f(x) > 0, so ln(f(x)) is defined
        samples = exp(samples - maxSample);

        for (int n = 0; n < numSample; ++n) {
            const T xy = samples[n];
            var.toNextVariance(mean, n, xy);

            const auto l = xy.value().squaredNorm();
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indices[n * getDim() + i];
                lossMat[index, i].toNextMean(counts[index, i], l);
                counts[index, i] += 1;
            }
        }
        lossMat.clamp_min(std::numeric_limits<Trv>::min());

        mean = ln(mean) + maxSample;
        var = ln(var + Trv(std::numeric_limits<T>::min()));
        var.value() = fma(Trv(2), maxSample, var.value());
    }
}
