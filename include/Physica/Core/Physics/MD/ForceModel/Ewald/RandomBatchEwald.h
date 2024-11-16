/*
 * Copyright 2023-2024 Weibo He.
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

#include "RSpaceEwald.h"

namespace Physica::Core {
    /**
     * TODO: Anisotropic condition
     * 
     * Reference:
     * [1] SIAM J. Sci. Comput. 43, B937-B960 (2021); https://doi.org/10.1137/20M1371385
     */
    template<Scalar T, class RandomType>
    class RandomBatchEwald : private RSpaceEwald<T, false> {
        using This = RandomBatchEwald<T, RandomType>;
        using Base = RSpaceEwald<T, false>;
        using Base::Dim;
        using typename Base::ValueType;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using SamplePool = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Dynamic, Dim>;
    public:
        using typename Base::BornChargeArray;
    private:
        T sumGauss;
        SamplePool samplePool;
        size_t batchSize;
    public:
        RandomBatchEwald() = default;
        RandomBatchEwald(size_t samplePoolSize, size_t batchSize_);
        RandomBatchEwald(size_t samplePoolSize, size_t batchSize_, LatticeMatrix lattice, VectorND<T> charges);
        RandomBatchEwald(const RandomBatchEwald&) = default;
        RandomBatchEwald(RandomBatchEwald&&) noexcept = default;
        ~RandomBatchEwald() = default;
        /* Operators */
        RandomBatchEwald& operator=(RandomBatchEwald obj) noexcept { swap(obj); return *this; }
        RandomBatchEwald& operator=(Base base);
        /* Operations */
        template<class Executor>
        [[nodiscard]] VectorND<T> force(const PositionMatrix& pos) const;
        using Base::force_short;
        template<class Executor>
        [[nodiscard]] VectorND<T> force_long(const PositionMatrix& pos) const;

        [[nodiscard]] T calcDefaultIntegralLimit() const;
        void swap(RandomBatchEwald& __restrict obj) noexcept;
        /* Getters */
        using Base::getNumParticle;
        using Base::getLattice;
        [[nodiscard]] size_t getSamplePoolSize() const noexcept { return samplePool.getRow(); }
        /* Setters */
        void setLattice(LatticeMatrix lattice);
        void setIntegralLimit(T integralLimit);
    private:
        void updateSumGauss();
        void updateSamplePool();
        inline void monteCarloStep(int& sample, T& prop_last, T deviation, T factor, T propAtZero) const;
        [[nodiscard]] Vector3D<T> randWaveG() const;
        [[nodiscard]] static bool checkParam(const LatticeMatrix& lattice);
    };

    template<Scalar T, class RandomType>
    RandomBatchEwald<T, RandomType>::RandomBatchEwald(size_t samplePoolSize, size_t batchSize_)
            : Base()
            , samplePool(samplePoolSize, Dim)
            , batchSize(batchSize_)  {}

    template<Scalar T, class RandomType>
    RandomBatchEwald<T, RandomType>::RandomBatchEwald(
            size_t samplePoolSize,
            size_t batchSize_,
            LatticeMatrix lattice,
            VectorND<T> charges) : This(samplePoolSize, batchSize_) {
        assert(checkParam(lattice) && "[Error]: Non-orthogonal lattice is not implemented");
        Base::setCharges(std::move(charges));
        Base::setLattice(std::move(lattice));
        setIntegralLimit(calcDefaultIntegralLimit());
    }

    template<Scalar T, class RandomType>
    RandomBatchEwald<T, RandomType>& RandomBatchEwald<T, RandomType>::operator=(Base base) {
        Base::swap(base);
        setIntegralLimit(calcDefaultIntegralLimit());
        return *this;
    }

    template<Scalar T, class RandomType>
    template<class Executor>
    VectorND<T> RandomBatchEwald<T, RandomType>::force(const PositionMatrix& pos) const {
        VectorND<T> result;
        auto kSpaceFuture = Executor::schedule([this, pos, &result]() {
            result = force_long<SequentialExecutor>(pos);
        });
        const VectorND<T> rSpaceSum = Base::template force_short<SequentialExecutor>(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }
    //TODO: time reversal symmetry?
    template<Scalar T, class RandomType>
    template<class Executor>
    VectorND<T> RandomBatchEwald<T, RandomType>::force_long(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        VectorND<T> kSpaceSum(numParticle * Dim, 0);
        VectorND<T> dots(numParticle);
        VectorND<T> sin_vec(numParticle);
        VectorND<T> cos_vec(numParticle);
        size_t numLoop = batchSize;
        for (size_t _ = 0; _ < numLoop; ++_) {
            const Vector3D<T> waveG = randWaveG();
            const T squaredNorm = T(waveG.squaredNorm());
            const bool isGammaPoint = squaredNorm < T(std::numeric_limits<T>::min());
            if (isGammaPoint) {
                numLoop += 1;
                continue;
            }
            dots = pos * waveG;
            sincos(dots, sin_vec, cos_vec);
            const auto& charges = Base::getCharges();
            const T sum_cos = cos_vec * charges;
            const T sum_sin = sin_vec * charges;
            const T factor = reciprocal(squaredNorm);
            for (size_t i = 0; i < numParticle; ++i) {
                auto force_i = kSpaceSum.template segment<3>(i * Dim, (i + 1) * Dim);
                const T temp = (sin_vec[i] * sum_cos - cos_vec[i] * sum_sin) * (factor * charges[i]);
                force_i[0] += temp * waveG[0];
                force_i[1] += temp * waveG[1];
                force_i[2] += temp * waveG[2];
            }
        }
        if constexpr (T::isReverseDiff)
            kSpaceSum.makeContinuous();
        kSpaceSum *= T(4 * M_PI) * Base::getInvVolume() * sumGauss / T(batchSize);
        return kSpaceSum;
    }

    template<Scalar T, class RandomType>
    T RandomBatchEwald<T, RandomType>::calcDefaultIntegralLimit() const {
        const T averageCellSize = cbrt(Base::getVolume());
        const T result = cbrt(T(getNumParticle())) * T(M_PI * 0.5) / averageCellSize; // Constant estimated from results of [1]
        return result;
    }

    template<Scalar T, class RandomType>
    void RandomBatchEwald<T, RandomType>::swap(RandomBatchEwald& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sumGauss.swap(obj.sumGauss);
        std::swap(batchSize, obj.batchSize);
        samplePool.swap(obj.samplePool);
    }

    template<Scalar T, class RandomType>
    void RandomBatchEwald<T, RandomType>::setLattice(LatticeMatrix lattice) {
        Base::setLattice(std::move(lattice));
        setIntegralLimit(calcDefaultIntegralLimit());
    }

    template<Scalar T, class RandomType>
    void RandomBatchEwald<T, RandomType>::setIntegralLimit(T integralLimit) {
        Base::setIntegralLimit(integralLimit);
        updateSumGauss();
        updateSamplePool();
    }

    template<Scalar T, class RandomType>
    void RandomBatchEwald<T, RandomType>::updateSumGauss() {
        sumGauss = 0;
        const T factor = reciprocal(square(ValueType(2) * Base::getIntegralLimit()));
        PeriodicCell<T, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            Base::getKSpaceSumRange(), Base::getRepLattice(), [this, factor](Vector3D<T> delta) {
                const T squaredNorm = T(delta.squaredNorm());
                const bool isNotGammaPoint = T(std::numeric_limits<T>::min()) < squaredNorm;
                if (isNotGammaPoint)
                    sumGauss += exp(-squaredNorm * factor);
            });
        sumGauss *= ValueType(2);
    }

    template<Scalar T, class RandomType>
    void RandomBatchEwald<T, RandomType>::updateSamplePool() {
        const size_t poolSize = getSamplePoolSize();
        for (size_t i = 0; i < Dim; ++i) {
            const auto& repLattice = Base::getRepLattice();
            const auto lattVec = repLattice.row(i);
            const T deviation = sqrt(T(2) * square(Base::getIntegralLimit() / lattVec.squaredNorm()));
            const T factor = reciprocal(deviation * M_SQRT2);
            const T propAtZero = erf(T(0.5) * factor);
            int sample = 0;
            T prop_last = reciprocal(propAtZero);
            for (size_t n = 0; n < poolSize; ++n) {
                monteCarloStep(sample, prop_last, deviation, factor, propAtZero);
                samplePool(n, i) = T(sample) * lattVec[i];
            }
        }
    }

    template<Scalar T, class RandomType>
    inline void RandomBatchEwald<T, RandomType>::monteCarloStep(
            int& sample,
            T& prop_last,
            T deviation,
            T factor,
            T propAtZero) const {
        auto& pool = RandomType::getInstance();
        const T rand = T::random_normal(pool);
        const T x0 = rand * deviation + T(rand.isPositive() ? 0.5 : -0.5);
        const int sample_new = double(x0);
        const T x1 = T(sample_new);

        T prop_new = exp(square(x1 / deviation) * T(-0.5));
        if (sample_new == 0)
            prop_new /= propAtZero;
        else {
            const T abs_x = abs(x1);
            prop_new /= (erf(factor * (abs_x + 0.5)) - erf(factor * (abs_x - 0.5))) * 0.5;
        }

        if (prop_new > prop_last) {
            sample = sample_new;
            prop_last = prop_new;
        }
        else {
            const T temp = T::random_uniform(pool);
            if (prop_new / prop_last > temp) {
                sample = sample_new;
                prop_last = prop_new;
            }
        }
    }

    template<Scalar T, class RandomType>
    Vector3D<T> RandomBatchEwald<T, RandomType>::randWaveG() const {
        auto dist = std::uniform_int_distribution<size_t>(0, getSamplePoolSize() - 1);
        auto& gen = RandomType::getInstance().getGen();
        Vector3D<T> result{};
        for (size_t i = 0; i < Dim; ++i)
            result[i] = samplePool(dist(gen), i);
        return result;
    }

    template<Scalar T, class RandomType>
    bool RandomBatchEwald<T, RandomType>::checkParam(const LatticeMatrix& lattice) {
        for (int i = 0; i < Dim; ++i) {
            for (int j = 0; j < Dim; ++j) {
                if (i == j)
                    continue;
                if (!lattice(i, j).isZero())
                    return false;
            }
        }
        return true;
    }
}

namespace Physica {
    using namespace Core;

    template<Scalar T, class RandomType>
    class Traits<RandomBatchEwald<T, RandomType>> : public Traits<RSpaceEwald<T, false>> {
    public:
        using REwaldType = RSpaceEwald<T, false>;
    };
}
