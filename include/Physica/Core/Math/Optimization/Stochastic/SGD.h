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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"
#include "Physica/Core/MultiPrecision/Differentiable.h"

namespace Physica::Core {
    /**
     * \class SGD Stochastic gradient descent using auto diff
     */
    template<class ScalarType>
    class SGD {
        static_assert(ScalarType::isDifferentiable, "[Error]: ScalarType must be differentiable");
        static_assert(!is_device_obj<ScalarType>::value, "[Error]: Include corresponding *.cuh file to enable CUDA support");
        using TracerType = typename ScalarType::TracerType;
    public:
        using PlainScalar = typename ScalarType::PlainScalar;
    private:
        constexpr static int AnyValue = 0;
    protected:
        PlainScalar learnRate;
        PlainScalar meanLearnRate;
        unsigned int batchSize;
        ScalarType from;
        ScalarType to;
    public:
        SGD(PlainScalar learnRate_, unsigned int batchSize_);
        SGD(const SGD&) = default;
        SGD(SGD&&) noexcept = default;
        ~SGD() = default;
        /* Operators */
        SGD& operator=(SGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void recordBegin();
        inline void recordEnd();
        void step() const;
        inline void zero_grad() const;
        void swap(SGD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] PlainScalar getLearnRate() const noexcept { return learnRate; }
        [[nodiscard]] PlainScalar getMeanLearnRate() const noexcept { return meanLearnRate; }
        [[nodiscard]] unsigned int getBatchSize() const noexcept { return batchSize; }
        /* Setters */
        void setLearnRate(PlainScalar lr);
    };

    template<class ScalarType>
    SGD<ScalarType>::SGD(PlainScalar learnRate_, unsigned int batchSize_) : batchSize(batchSize_) {
        assert(batchSize > 0);
        setLearnRate(learnRate_);
    }

    template<class ScalarType>
    inline void SGD<ScalarType>::recordBegin() {
        to = ScalarType(AnyValue);
    }
    
    template<class ScalarType>
    inline void SGD<ScalarType>::recordEnd() {
        from = ScalarType(AnyValue);
    }

    template<class ScalarType>
    void SGD<ScalarType>::step() const {
        using SegmentType = typename TracerType::SegmentType;
        using DiffScalar = typename TracerType::DiffScalar;
        auto& tracer = TracerType::getInstance();
        tracer.forSegmentInRange(from, to, [this](SegmentType& segment) {
            segment.forNodeInRange(from, to, [this](DiffScalar s) {
                s.setValue(s.getValue() - meanLearnRate * s.getGrad());
            });
        });
    }

    template<class ScalarType>
    inline void SGD<ScalarType>::zero_grad() const {
        TracerType::getInstance().zero_grad(from, to);
    }

    template<class ScalarType>
    void SGD<ScalarType>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
        meanLearnRate.swap(obj.meanLearnRate);
        std::swap(batchSize, obj.batchSize);
        from.swap(obj.from);
        to.swap(obj.to);
    }

    template<class PlainScalar>
    void SGD<PlainScalar>::setLearnRate(PlainScalar lr) {
        assert(!lr.isZero() && "[Error]: 0 learn rate does nothing");
        learnRate = lr;
        meanLearnRate = lr / PlainScalar(batchSize);
    }
}
