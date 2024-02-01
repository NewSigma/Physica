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

#include <memory>
#include "SGD.h"

namespace Physica::Core {
    template<class ScalarType>
    class MomentumSGD : public SGD<ScalarType> {
        static_assert(ScalarType::isDifferentiable, "[Error]: ScalarType must be differentiable");
        static_assert(!Utils::is_device_obj<ScalarType>::value, "[Error]: Not implemented");
        using Base = SGD<ScalarType>;
        using typename Base::PlainScalar;
        using TracerType = DiffTracer<PlainScalar>;
    private:
        using Base::from;
        using Base::to;
        Vector<PlainScalar> lastGrad;
        PlainScalar momentum;
    public:
        MomentumSGD(PlainScalar momentum_, PlainScalar learnRate, unsigned int batchSize);
        MomentumSGD(const MomentumSGD&) = default;
        MomentumSGD(MomentumSGD&&) noexcept = default;
        ~MomentumSGD() = default;
        /* Operators */
        MomentumSGD& operator=(MomentumSGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void recordEnd();
        void step();
        void swap(MomentumSGD& __restrict obj) noexcept;
    private:
        using Base::step;
    };

    template<class ScalarType>
    MomentumSGD<ScalarType>::MomentumSGD(PlainScalar momentum_, PlainScalar learnRate, unsigned int batchSize)
            : Base(std::move(learnRate), batchSize)
            , momentum(std::move(momentum_)) {
        assert(momentum.isPositive() && "[Error]: Invalid momentum");
    }

    template<class ScalarType>
    void MomentumSGD<ScalarType>::recordEnd() {
        Base::recordEnd();
        lastGrad.resize(TracerType::distance(from, to), PlainScalar(0));
    }

    template<class ScalarType>
    void MomentumSGD<ScalarType>::step() {
        using SegmentType = typename TracerType::SegmentType;
        using DiffScalar = typename TracerType::DiffScalar;
        auto& tracer = TracerType::getInstance();
        size_t i = 0;
        tracer.forSegmentInRange(from, to, [this, &i](SegmentType& segment) {
            segment.forNodeInRange(from, to, [this, &i](DiffScalar s) {
                lastGrad[i] = momentum * lastGrad[i] + s.getGrad();
                s.setValue(s.getValue() - Base::getMeanLearnRate() * lastGrad[i]);
                i += 1;
            });
        });
    }

    template<class ScalarType>
    void MomentumSGD<ScalarType>::swap(MomentumSGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        momentum.swap(obj.momentum);
    }
}
