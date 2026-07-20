/*
 * Copyright 2025-2026 Weibo He.
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
    Vegas<T, TakeLn>::Vegas(Cube cube, int numSample, int numPoint, Trv compressRate, Trv lr, Trv momentum)
            : cube(std::move(cube))
            , pointGrid(numPoint, getDim())
            , compressRate(compressRate)
            , lr(lr)
            , momentum(momentum)
            , losses(numPoint - 1, getDim())
            , counts(numPoint - 1, getDim())
            , indices(getDim(), numSample) {
        assert(numPoint > 2 && "[Error]: Invalid point number");
        assert(compressRate.isPositive() && "[Error]: Rate = 0 implies no grid refinement");
        assert(compressRate < Trv(2) && "[Error]: Rate should be ~ 1");
        assert(lr.isPositive() && (lr <= Trv(1)) && "[Error]: Invalid param");
        assert(!momentum.isNegative() && (momentum <= Trv(1)) && "[Error]: Invalid param");
        for (size_t i = 0; i < getDim(); ++i)
            mesh_uniform(i);
    }

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::mesh_uniform(size_t dim) {
        pointGrid.col(dim) = VectorND<Trv>::linspace(cube[0, dim], cube[1, dim], getNumPoint());
    }

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::mesh_tanh(size_t dim, Trv range) {
        assert(range.isPositive());
        const auto k = (cube[1, dim] - cube[0, dim]) * 0.5;
        const auto b = (cube[1, dim] + cube[0, dim]) * 0.5;
        auto p = VectorND<Trv>::linspace(-range, range, getNumPoint());
        p = tanh(p);
        p[0] = -1;
        p[getNumPoint() - 1] = 1;
        pointGrid.col(dim) = k * p + b;
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    auto Vegas<T, TakeLn>::warmup(std::invocable<VectorND<Trv>> auto fn, int numWarm) -> Trv {
        assert(numWarm >= 0 && "[Error]: Invalid param");
        checkIntegrand<VectorND<Trv>>(fn);

        DenseMatrix<Trv> losses_old;
        bool hasMomentum = momentum.isPositive();
        if (hasMomentum) {
            losses_old.resize(losses);
            losses_old.zeros();
        }

        for (int i = 0; i < numWarm; ++i) {
            init<R>();
            if constexpr (TakeLn)
                sample_ln<R, P>(fn);
            else
                sample<R, P>(fn);

            compress<P>();
            if (hasMomentum) {
                if (i != 0)
                    losses = losses * (Trv(1) - momentum) + losses_old * momentum;
                refineGrid<P>();
                losses_old = losses;
            }
            else
                refineGrid<P>();
        }
        return calcGridLossImpl();
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    AdaptiveProcess<T, TakeLn> Vegas<T, TakeLn>::integral(std::invocable<VectorND<Trv>> auto fn, int numRefine) {
        checkIntegrand<VectorND<Trv>>(fn);

        VectorND<T> means(numRefine);
        VectorND<T> vars(numRefine);
        VectorND<Trv> loss(numRefine);
        DenseMatrix<Trv> losses_old;
        bool hasMomentum = momentum.isPositive();
        if (hasMomentum) {
            losses_old.resize(losses);
            losses_old.zeros();
        }

        for (int i = 0; i < numRefine; ++i) {
            init<R>();
            if constexpr (TakeLn) {
                auto [mean, var] = sample_ln<R, P>(fn);
                means[i] = mean;
                vars[i] = var;
            }
            else {
                auto [mean, var] = sample<R, P>(fn);
                means[i] = mean;
                vars[i] = var;
            }
            loss[i] = calcGridLossImpl();

            compress<P>();
            if (hasMomentum) {
                if (i != 0)
                    losses = losses * (Trv(1) - momentum) + losses_old * momentum;
                refineGrid<P>();
                losses_old = losses;
            }
            else
                refineGrid<P>();
        }
        return AdaptiveProcess<T, TakeLn>(std::move(means), std::move(vars), std::move(loss), getNumSample());
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    auto Vegas<T, TakeLn>::calcGridLoss(std::invocable<VectorND<Trv>> auto fn) const -> Trv {
        init<R>();
        if constexpr (TakeLn)
            sample_ln<R, P>(fn);
        else
            sample<R, P>(fn);
        return calcGridLossImpl();
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, bool TakeLn>
    const H5Group Vegas<T, TakeLn>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        group.readAttr("CompressRate", compressRate);
        group.readAttr("lr", lr);

        cube.read(group, "Cube");
        pointGrid.read(group, "Grid");
        losses.resize(getNumPoint() - 1, cube.getCol());
        counts.resize(getNumPoint() - 1, cube.getCol());
        return group;
    }

    template<Scalar T, bool TakeLn>
    H5Group Vegas<T, TakeLn>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        group.writeAttr("CompressRate", compressRate);
        group.writeAttr("lr", lr);

        cube.write(group, "Cube");
        pointGrid.write(group, "Grid");
        return group;
    }
#endif

    template<Scalar T, bool TakeLn>
    void Vegas<T, TakeLn>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cube.swap(obj.cube);
        pointGrid.swap(obj.pointGrid);
        compressRate.swap(obj.compressRate);
        lr.swap(obj.lr);
        momentum.swap(obj.momentum);
        losses.swap(obj.losses);
        counts.swap(obj.counts);
    }

    template<Scalar T, bool TakeLn>
    template<ExecutePolicy P>
    void Vegas<T, TakeLn>::compress() {
        parallel_for<P>([this](size_t dim) {
            auto loss = losses.col(dim);
            // Normalized values fall into a range that is feasible for compression function.
            const VectorND<Trv> buffer = loss * reciprocal(loss.sum());
            const Vector3D<Trv> kernel{1.0 / 8, 6.0 / 8, 1.0 / 8};
            size_t i = 0;
            loss[0] = fma(Trv(7.0 / 8), buffer[0], Trv(1.0 / 8) * buffer[1]);
            for (; i < loss.getLength() - 2; ++i)
                loss[i + 1] = buffer.template segment<3>(i, i + 3) * kernel;
            loss[i + 1] = fma(Trv(7.0 / 8), buffer[i + 1], Trv(1.0 / 8) * buffer[i]);

            loss = pow(divide(loss - Trv(1), ln(loss + Trv(std::numeric_limits<T>::min()))), compressRate);
        }, getDim(), 0).wait();
    }

    template<Scalar T, bool TakeLn>
    template<RNG R>
    void Vegas<T, TakeLn>::init() {
        // Loss has minimal value to avoid the grid size reducing to 0
        losses = Trv(std::numeric_limits<T>::min());
        counts.zeros();
        R::getInstance().random_int(indices.asArray(), 0, getNumPoint() - 2);
    }

    template<Scalar T, bool TakeLn>
    template<ExecutePolicy P>
    void Vegas<T, TakeLn>::refineGrid() {
        parallel_for<P>([this](size_t dim) {
            const auto mean = losses.col(dim).mean();
            auto oldPoints = pointGrid.col(dim);
            VectorND<Trv> newPoints(getNumPoint());
            newPoints.front() = oldPoints.front();
            Trv temp = 0;
            for (size_t i = 1, j = 0; i < newPoints.getLength() - 1; ++i) {
                while (temp < mean) {
                    assert(j < getNumPoint() && "[Error]: Unexpected not enough loss, this is likely a bug");
                    temp += losses[j, dim];
                    j += 1;
                }
                temp -= mean;
                Trv delta = oldPoints[j] - oldPoints[j - 1];
                newPoints[i] = fma(temp / losses[j - 1, dim], -delta, oldPoints[j]);
            }
            newPoints.back() = oldPoints.back();
            oldPoints = newPoints * lr + oldPoints * (Trv(1) - lr);
        }, getDim(), 0).wait();
    }

    template<Scalar T, bool TakeLn>
    auto Vegas<T, TakeLn>::calcGridLossImpl() const -> Trv {
        Trv sumvar = 0;
        for (size_t i = 0; i < getDim(); ++i) {
            const auto col = losses.col(i);
            const Trv prior = col.mean();
            sumvar += sqrt(col.variance(prior) / prior.squaredNorm() / Trv(getNumPoint()));
        }
        return sumvar / Trv(getDim());
    }

    template<Scalar T, bool TakeLn>
    template<RNG R>
    auto Vegas<T, TakeLn>::transform(std::span<int> indices) const -> Array<VectorND<Trv>, 2> {
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
        return {std::move(x), std::move(deltas)};
    }

    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    Vector2D<T> Vegas<T, TakeLn>::sample(std::invocable<VectorND<Trv>> auto fn) {
        const size_t numSample = getNumSample();
        VectorND<T> samples(numSample);
        parallel_for<P>([&, this](size_t n) {
            const auto [x, deltas] = transform<R>(indices.col(n));
            const T y = fn(x);
            assert(y.isFinite() && "[Error]: Bad value");

            const T xy = y * (deltas * Trv(getNumPoint())).prod();
            samples[n] = xy;
        }, numSample, 0).wait();
        return statistic(samples);
    }
    /**
     * Assuming ln(f(x)) is defined
     */
    template<Scalar T, bool TakeLn>
    template<RNG R, ExecutePolicy P>
    Vector2D<T> Vegas<T, TakeLn>::sample_ln(std::invocable<VectorND<Trv>> auto fn) {
        const size_t numSample = getNumSample();
        VectorND<T> samples(numSample);
        parallel_for<P>([&, this](size_t n) {
            const auto [x, deltas] = transform<R>(indices.col(n));
            const T lny = fn(x);
            assert(lny.isFinite() && "[Error]: Bad value");

            T lnxy = lny + ln(deltas).sum();
            lnxy.value().real() = fma(Trv(getDim()), ln(Trv(getNumPoint())), lnxy.value().real());
            samples[n] = lnxy;
        }, numSample, 0).wait();

        const Trv maximum = [&]() {
            if constexpr (T::isComplex())
                return samples.reals().max().value();
            else
                return samples.max().value();
        }();
        samples = exp(samples - maximum);

        auto [mean, var] = statistic(samples);
        mean = ln(mean) + maximum;
        var = ln(var + Trv(std::numeric_limits<T>::min()));
        var.value().real() = fma(Trv(2), maximum, var.value().real());
        assert(var.isFinite() && "[Error]: maximum is subnormal?");
        return {mean, var};
    }

    template<Scalar T, bool TakeLn>
    Vector2D<T> Vegas<T, TakeLn>::statistic(const VectorND<T>& samples) {
        T mean, var;
        for (size_t n = 0; n < samples.getLength(); ++n) {
            const T xy = samples[n];
            var.toNextVariance(mean, n, xy);

            const auto l = xy.value().squaredNorm();
            for (size_t i = 0; i < getDim(); ++i) {
                const auto index = indices[i, n];
                losses[index, i].toNextMean(counts[index, i], l);
                counts[index, i] += 1;
            }
        }
        losses.clamp_min(std::numeric_limits<Trv>::min());
        return {mean, var};
    }

    template<Scalar T, bool TakeLn>
    template<class Sample>
    consteval void Vegas<T, TakeLn>::checkIntegrand(std::invocable<Sample> auto& fn) noexcept {
        using R = std::invoke_result<decltype(fn), Sample>::type;
        static_assert(std::is_same<R, T>::value, "[Error]: Return type does not match Vegas's working type");
    }
}
