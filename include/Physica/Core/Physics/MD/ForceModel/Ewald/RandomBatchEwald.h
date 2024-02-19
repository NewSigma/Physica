/*
 * Copyright 2023-2024 WeiBo He.
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
    template<class ScalarType, class RandomPoolType> class RandomBatchEwald;

    namespace Internal {
        template<class ScalarType, class RandomPoolType>
        class Traits<RandomBatchEwald<ScalarType, RandomPoolType>> : public Traits<RSpaceEwald<ScalarType, false>> {
        public:
            using REwaldType = RSpaceEwald<ScalarType, false>;
        };
    }
    /**
     * TODO: Anisotropic condition
     * 
     * Reference:
     * [1] SIAM J. Sci. Comput. 43, B937-B960 (2021); https://doi.org/10.1137/20M1371385
     */
    template<class ScalarType, class RandomPoolType>
    class RandomBatchEwald : private RSpaceEwald<ScalarType, false> {
        using This = RandomBatchEwald<ScalarType, RandomPoolType>;
        using Base = RSpaceEwald<ScalarType, false>;
        using Base::Dim;
        using typename Base::PlainScalar;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Vector3D;
        using SamplePool = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Element, Dynamic, Dim>;
    public:
        using typename Base::BornChargeArray;
    private:
        ScalarType sumGauss;
        SamplePool samplePool;
        size_t batchSize;
    public:
        RandomBatchEwald() = default;
        RandomBatchEwald(size_t samplePoolSize, size_t batchSize_);
        RandomBatchEwald(size_t samplePoolSize, size_t batchSize_, LatticeMatrix lattice, Vector<ScalarType> charges);
        RandomBatchEwald(const RandomBatchEwald&) = default;
        RandomBatchEwald(RandomBatchEwald&&) noexcept = default;
        ~RandomBatchEwald() = default;
        /* Operators */
        RandomBatchEwald& operator=(RandomBatchEwald obj) noexcept { swap(obj); return *this; }
        RandomBatchEwald& operator=(Base base);
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const PositionMatrix& pos) const;
        using Base::force_short;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const PositionMatrix& pos) const;

        [[nodiscard]] ScalarType calcDefaultIntegralLimit() const;
        void swap(RandomBatchEwald& __restrict obj) noexcept;
        /* Getters */
        using Base::getNumParticle;
        using Base::getLattice;
        [[nodiscard]] size_t getSamplePoolSize() const noexcept { return samplePool.getRow(); }
        /* Setters */
        void setLattice(LatticeMatrix lattice);
        void setIntegralLimit(ScalarType integralLimit);
    private:
        void updateSumGauss();
        void updateSamplePool();
        inline void monteCarloStep(int& sample, ScalarType& prop_last, ScalarType deviation, ScalarType factor, ScalarType propAtZero) const;
        [[nodiscard]] Vector3D randWaveG() const;
        [[nodiscard]] static bool checkParam(const LatticeMatrix& lattice);
    };

    template<class ScalarType, class RandomPoolType>
    RandomBatchEwald<ScalarType, RandomPoolType>::RandomBatchEwald(size_t samplePoolSize, size_t batchSize_)
            : Base()
            , samplePool(samplePoolSize, Dim)
            , batchSize(batchSize_)  {}

    template<class ScalarType, class RandomPoolType>
    RandomBatchEwald<ScalarType, RandomPoolType>::RandomBatchEwald(
            size_t samplePoolSize,
            size_t batchSize_,
            LatticeMatrix lattice,
            Vector<ScalarType> charges) : This(samplePoolSize, batchSize_) {
        assert(checkParam(lattice) && "[Error]: Non-orthogonal lattice is not implemented");
        Base::setCharges(std::move(charges));
        Base::setLattice(std::move(lattice));
        setIntegralLimit(calcDefaultIntegralLimit());
    }

    template<class ScalarType, class RandomPoolType>
    RandomBatchEwald<ScalarType, RandomPoolType>& RandomBatchEwald<ScalarType, RandomPoolType>::operator=(Base base) {
        Base::swap(base);
        setIntegralLimit(calcDefaultIntegralLimit());
        return *this;
    }

    template<class ScalarType, class RandomPoolType>
    template<class Executor>
    Vector<ScalarType> RandomBatchEwald<ScalarType, RandomPoolType>::force(const PositionMatrix& pos) const {
        Vector<ScalarType> result;
        auto kSpaceFuture = Executor::schedule([this, pos, &result]() {
            result = force_long<SequentialExecutor>(pos);
        });
        const Vector<ScalarType> rSpaceSum = Base::template force_short<SequentialExecutor>(pos);
        Executor::auto_wait(kSpaceFuture);
        result += rSpaceSum;
        return result;
    }
    //TODO: time reversal symmetry?
    template<class ScalarType, class RandomPoolType>
    template<class Executor>
    Vector<ScalarType> RandomBatchEwald<ScalarType, RandomPoolType>::force_long(const PositionMatrix& pos) const {
        const size_t numParticle = getNumParticle();
        Vector<ScalarType> kSpaceSum(numParticle * Dim, 0);
        Vector<ScalarType> dots(numParticle);
        Vector<ScalarType> sin_vec(numParticle);
        Vector<ScalarType> cos_vec(numParticle);
        size_t numLoop = batchSize;
        for (size_t _ = 0; _ < numLoop; ++_) {
            const Vector3D waveG = randWaveG();
            const ScalarType squaredNorm = ScalarType(waveG.squaredNorm());
            const bool isGammaPoint = squaredNorm < ScalarType(std::numeric_limits<ScalarType>::min());
            if (isGammaPoint) {
                numLoop += 1;
                continue;
            }
            dots = pos * waveG;
            sincos(dots, sin_vec, cos_vec);
            const auto& charges = Base::getCharges();
            const ScalarType sum_cos = cos_vec * charges;
            const ScalarType sum_sin = sin_vec * charges;
            const ScalarType factor = reciprocal(squaredNorm);
            for (size_t i = 0; i < numParticle; ++i) {
                auto force_i = kSpaceSum.template segment<3>(i * Dim, (i + 1) * Dim);
                const ScalarType temp = (sin_vec[i] * sum_cos - cos_vec[i] * sum_sin) * (factor * charges[i]);
                force_i[0] += temp * waveG[0];
                force_i[1] += temp * waveG[1];
                force_i[2] += temp * waveG[2];
            }
        }
        if constexpr (ScalarType::isReverseDiff)
            kSpaceSum.makeContinuous();
        kSpaceSum *= ScalarType(4 * M_PI) * Base::getInvVolume() * sumGauss / ScalarType(batchSize);
        return kSpaceSum;
    }

    template<class ScalarType, class RandomPoolType>
    ScalarType RandomBatchEwald<ScalarType, RandomPoolType>::calcDefaultIntegralLimit() const {
        const ScalarType averageCellSize = cbrt(Base::getVolume());
        const ScalarType result = cbrt(ScalarType(getNumParticle())) * ScalarType(M_PI * 0.5) / averageCellSize; // Constant estimated from results of [1]
        return result;
    }

    template<class ScalarType, class RandomPoolType>
    void RandomBatchEwald<ScalarType, RandomPoolType>::swap(RandomBatchEwald& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        sumGauss.swap(obj.sumGauss);
        std::swap(batchSize, obj.batchSize);
        samplePool.swap(obj.samplePool);
    }

    template<class ScalarType, class RandomPoolType>
    void RandomBatchEwald<ScalarType, RandomPoolType>::setLattice(LatticeMatrix lattice) {
        Base::setLattice(std::move(lattice));
        setIntegralLimit(calcDefaultIntegralLimit());
    }

    template<class ScalarType, class RandomPoolType>
    void RandomBatchEwald<ScalarType, RandomPoolType>::setIntegralLimit(ScalarType integralLimit) {
        Base::setIntegralLimit(integralLimit);
        updateSumGauss();
        updateSamplePool();
    }

    template<class ScalarType, class RandomPoolType>
    void RandomBatchEwald<ScalarType, RandomPoolType>::updateSumGauss() {
        sumGauss = 0;
        const ScalarType factor = reciprocal(square(PlainScalar(2) * Base::getIntegralLimit()));
        PeriodicCell<ScalarType, Dim>::forReducedCellInRange( // Reduce cell using time reversal symmetry
            Base::getKSpaceSumRange(), Base::getRepLattice(), [this, factor](Vector3D delta) {
                const ScalarType squaredNorm = ScalarType(delta.squaredNorm());
                const bool isNotGammaPoint = ScalarType(std::numeric_limits<ScalarType>::min()) < squaredNorm;
                if (isNotGammaPoint)
                    sumGauss += exp(-squaredNorm * factor);
            });
        sumGauss *= PlainScalar(2);
    }

    template<class ScalarType, class RandomPoolType>
    void RandomBatchEwald<ScalarType, RandomPoolType>::updateSamplePool() {
        const size_t poolSize = getSamplePoolSize();
        for (size_t i = 0; i < Dim; ++i) {
            const auto& repLattice = Base::getRepLattice();
            const auto lattVec = repLattice.row(i);
            const ScalarType deviation = sqrt(ScalarType(2) * square(Base::getIntegralLimit() / lattVec.squaredNorm()));
            const ScalarType factor = reciprocal(deviation * M_SQRT2);
            const ScalarType propAtZero = erf(ScalarType(0.5) * factor);
            int sample = 0;
            ScalarType prop_last = reciprocal(propAtZero);
            for (size_t n = 0; n < poolSize; ++n) {
                monteCarloStep(sample, prop_last, deviation, factor, propAtZero);
                samplePool(n, i) = ScalarType(sample) * lattVec[i];
            }
        }
    }

    template<class ScalarType, class RandomPoolType>
    inline void RandomBatchEwald<ScalarType, RandomPoolType>::monteCarloStep(
            int& sample,
            ScalarType& prop_last,
            ScalarType deviation,
            ScalarType factor,
            ScalarType propAtZero) const {
        auto& pool = RandomPoolType::getInstance();
        const ScalarType rand = ScalarType::random_normal(pool);
        const ScalarType x0 = rand * deviation + ScalarType(rand.isPositive() ? 0.5 : -0.5);
        const int sample_new = double(x0);
        const ScalarType x1 = ScalarType(sample_new);

        ScalarType prop_new = exp(square(x1 / deviation) * ScalarType(-0.5));
        if (sample_new == 0)
            prop_new /= propAtZero;
        else {
            const ScalarType abs_x = abs(x1);
            prop_new /= (erf(factor * (abs_x + 0.5)) - erf(factor * (abs_x - 0.5))) * 0.5;
        }

        if (prop_new > prop_last) {
            sample = sample_new;
            prop_last = prop_new;
        }
        else {
            const ScalarType temp = ScalarType::random_uniform(pool);
            if (prop_new / prop_last > temp) {
                sample = sample_new;
                prop_last = prop_new;
            }
        }
    }

    template<class ScalarType, class RandomPoolType>
    typename RandomBatchEwald<ScalarType, RandomPoolType>::Vector3D RandomBatchEwald<ScalarType, RandomPoolType>::randWaveG() const {
        auto dist = std::uniform_int_distribution<size_t>(0, getSamplePoolSize() - 1);
        auto& gen = RandomPoolType::getInstance().getGen();
        Vector3D result{};
        for (size_t i = 0; i < Dim; ++i)
            result[i] = samplePool(dist(gen), i);
        return result;
    }

    template<class ScalarType, class RandomPoolType>
    bool RandomBatchEwald<ScalarType, RandomPoolType>::checkParam(const LatticeMatrix& lattice) {
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
