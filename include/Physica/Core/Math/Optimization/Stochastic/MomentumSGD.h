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

#include <memory>
#include "SGD.h"

namespace Physica::Core {
    template<Scalar T>
    class MomentumSGD : public SGD<T> {
        static_assert(T::isDifferentiable, "[Error]: T must be differentiable");
        static_assert(!is_device_obj<T>::value, "[Error]: Not implemented");
        using Base = SGD<T>;
        using typename Base::ValueType;
        using TracerType = T::TracerType;
    private:
        using Base::from;
        using Base::to;
        VectorND<ValueType> lastGrad;
        ValueType momentum;
    public:
        MomentumSGD(ValueType momentum_, ValueType learnRate, unsigned int batchSize);
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

    template<Scalar T>
    MomentumSGD<T>::MomentumSGD(ValueType momentum_, ValueType learnRate, unsigned int batchSize)
            : Base(std::move(learnRate), batchSize)
            , momentum(std::move(momentum_)) {
        assert(momentum.isPositive() && "[Error]: Invalid momentum");
    }

    template<Scalar T>
    void MomentumSGD<T>::recordEnd() {
        Base::recordEnd();
        lastGrad.resize(TracerType::distance(from, to), ValueType(0));
    }

    template<Scalar T>
    void MomentumSGD<T>::step() {
        using SegmentType = TracerType::SegmentType;
        using DiffScalar = TracerType::DiffScalar;
        auto& tracer = TracerType::getInstance();
        size_t i = 0;
        tracer.forSegmentInRange(from, to, [this, &i](SegmentType& segment) {
            segment.forNodeInRange(from, to, [this, &i](DiffScalar s) {
                lastGrad[i] = momentum * lastGrad[i] + s.template getGrad<>();
                s.setValue(s.getValue() - Base::getMeanLearnRate() * lastGrad[i]);
                i += 1;
            });
        });
    }

    template<Scalar T>
    void MomentumSGD<T>::swap(MomentumSGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        momentum.swap(obj.momentum);
    }
}
