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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h>
#include "Hamilton/HamiltonMatrix.h"

namespace Physica::Core {
    /**
     * NVE ensemble referenced from [1]
     * NVT ensemble referenced from [2]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 108, 240401 (2012); https://doi.org/10.1103/PhysRevLett.108.240401
     * [2] Phys. Rev. Lett. 111, 010401 (2013); https://doi.org/10.1103/PhysRevLett.111.010401
     */
    template<class ScalarType>
    class TPQ : public VectorND<ScalarType> {
        using This = TPQ<ScalarType>;
        using Base = VectorND<ScalarType>;
        using RealType = typename ScalarType::RealType;

        RealType beta;
        RealType traceMu;
        int numMinCostTerm;
        int numSplit;
    public:
        TPQ() = default;
        TPQ(size_t length);
        TPQ(size_t length, RealType traceMu_, int numMinCostTerm_, int numSplit_);
        TPQ(const This&) = default;
        TPQ(This&&) noexcept = default;
        ~TPQ() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class ModelType, class Executor = SequentialExecutor>
        void pre_nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, RealType deltaBeta);
        template<class ModelType, class Executor = SequentialExecutor>
        void nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, RealType deltaBeta);

        [[nodiscard]] inline RealType calcPartitionXi() const;
        [[nodiscard]] inline RealType lnPartitionXi() const;
        [[nodiscard]] RealType lnSquaredDot(const VectorND<RealType>& other) const;
        void swap(This& __restrict obj) noexcept;

        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen, RealType norm = 0);
        /* Getters */
        [[nodiscard]] RealType getBeta() const noexcept { return beta; }
        [[nodiscard]] RealType getTraceMu() const noexcept { return traceMu; }
        [[nodiscard]] int getNumMinCostTerm() const noexcept { return numMinCostTerm; }
        [[nodiscard]] int getNumSplit() const noexcept { return numSplit; }
        [[nodiscard]] const Base& asVector() const noexcept { return *this; }
        [[nodiscard]] Base& asVector() noexcept { return *this; }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] static This random_normal(size_t len, RandomGenerator& gen, RealType norm = 0);
    private:
        using Base::random_uniform;
        using Base::random_any;
        [[nodiscard]] bool isPrepared() const;
    };

    template<class ScalarType>
    TPQ<ScalarType>::TPQ(size_t length)
            : Base(length)
            , beta(0)
            , traceMu(std::numeric_limits<RealType>::max())
            , numMinCostTerm(0)
            , numSplit(0) {
        assert(length > 0 && "[Error]: Invalid length");
    }

    template<class ScalarType>
    TPQ<ScalarType>::TPQ(size_t length, RealType traceMu_, int numMinCostTerm_, int numSplit_)
            : TPQ(length) {
        traceMu = std::move(traceMu_);
        numMinCostTerm = numMinCostTerm_;
        numSplit = numSplit_;
        assert(isPrepared() && "[Error]: Invalid params");
    }

    template<class ScalarType>
    template<class ModelType, class Executor>
    void TPQ<ScalarType>::pre_nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, RealType deltaBeta) {
        const RealType factor = deltaBeta * -0.5;
        const auto& hamiltonH = hamiltonH_.getDerived();
        const auto expr1 = factor * hamiltonH;
        const auto expr2 = exp(expr1);
        const auto expr3 = expr2 * (*this);
        traceMu = expr3.calcTraceMu();
        const auto params = expr3.template calcParam<Executor>(traceMu);
        numMinCostTerm = params.first;
        numSplit = params.second;
    }

    template<class ScalarType>
    template<class ModelType, class Executor>
    void TPQ<ScalarType>::nvt_step(const HamiltonMatrix<ModelType>& hamiltonH_, RealType deltaBeta) {
        if (deltaBeta.isZero())
            return;
        using BufferType = VectorND<ScalarType>;
        const RealType factor = deltaBeta * -0.5;
        const auto& hamiltonH = hamiltonH_.getDerived();
        const auto expr1 = factor * hamiltonH;
        const auto expr2 = exp(expr1);
        const auto expr3 = expr2 * (*this);
        BufferType dot(Base::getLength());
        if (isPrepared())
            expr3.template assignTo<BufferType, Executor>(dot, traceMu, std::make_pair(numMinCostTerm, numSplit));
        else
            expr3.template assignTo<BufferType, Executor>(dot);
        Base::swap(dot);
        beta += deltaBeta;
    }

    template<class ScalarType>
    inline typename TPQ<ScalarType>::RealType TPQ<ScalarType>::calcPartitionXi() const {
        return Base::squaredNorm();
    }

    template<class ScalarType>
    inline typename TPQ<ScalarType>::RealType TPQ<ScalarType>::lnPartitionXi() const {
        const bool isUnderflow = abs(*this).max() < RealType(std::numeric_limits<RealType>::min());
        if (isUnderflow)
            return RealType(-std::numeric_limits<ScalarType>::max());
        return Base::lnSquaredNorm();
    }

    template<class ScalarType>
    typename TPQ<ScalarType>::RealType TPQ<ScalarType>::lnSquaredDot(const VectorND<RealType>& other) const {
        assert(Base::getLength() == other.getLength() && "[Error]: Dimensions do not match");
        const RealType maxabs = abs(*this).max();
        const bool isUnderflow = maxabs < RealType(std::numeric_limits<RealType>::min());
        if (isUnderflow)
            return RealType(-std::numeric_limits<ScalarType>::max());
        const RealType factor = reciprocal(maxabs);
        const auto expr1 = (*this) * factor;
        const auto expr2 = toSquaredNormVector(expr1);
        const RealType dot = hadamard(expr2, other).sum();
        if (dot.isZero())
            return RealType(-std::numeric_limits<ScalarType>::max());
        return ln(dot) + RealType(2) * ln(maxabs);
    }

    template<class ScalarType>
    void TPQ<ScalarType>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        beta.swap(obj.beta);
        traceMu.swap(obj.traceMu);
        std::swap(numMinCostTerm, obj.numMinCostTerm);
        std::swap(numSplit, obj.numSplit);
    }

    template<class ScalarType>
    template<class RandomGenerator>
    inline void TPQ<ScalarType>::random_normal(RandomGenerator& gen, RealType norm) {
        const bool useDefault = norm.isZero();
        if (useDefault) // Default norm is as small as possible while keep all effective digits
            norm = RealType(std::numeric_limits<RealType>::min() / std::numeric_limits<RealType>::epsilon());
        Base::random_normal(gen);
        (*this) *= norm;
    }

    template<class ScalarType>
    template<class RandomGenerator>
    TPQ<ScalarType> TPQ<ScalarType>::random_normal(size_t len, RandomGenerator& gen, RealType norm) {
        This result(len);
        result.random_normal(gen, norm);
        return result;
    }

    template<class ScalarType>
    bool TPQ<ScalarType>::isPrepared() const {
        const bool muReady = traceMu != RealType(std::numeric_limits<RealType>::max());
        return muReady && numMinCostTerm > 0 && numSplit > 0;
    }
}
