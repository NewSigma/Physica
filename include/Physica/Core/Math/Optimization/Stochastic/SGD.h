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
#include "Physica/Core/MultiPrecision/Diff.h"

namespace Physica::Core {
    /**
     * \class SGD Stochastic gradient descent using auto diff
     */
    template<Scalar T>
    class SGD {
        static_assert(T::isDifferentiable, "[Error]: T must be differentiable");
        static_assert(!is_device_obj<T>::value, "[Error]: Include corresponding *.cuh file to enable CUDA support");
        using TracerType = typename T::TracerType;
    public:
        using ValueType = typename T::ValueType;
    private:
        constexpr static int AnyValue = 0;
    protected:
        ValueType learnRate;
        ValueType meanLearnRate;
        unsigned int batchSize;
        T from;
        T to;
    public:
        SGD(ValueType learnRate_, unsigned int batchSize_);
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
        [[nodiscard]] ValueType getLearnRate() const noexcept { return learnRate; }
        [[nodiscard]] ValueType getMeanLearnRate() const noexcept { return meanLearnRate; }
        [[nodiscard]] unsigned int getBatchSize() const noexcept { return batchSize; }
        /* Setters */
        void setLearnRate(ValueType lr);
    };

    template<Scalar T>
    SGD<T>::SGD(ValueType learnRate_, unsigned int batchSize_) : batchSize(batchSize_) {
        assert(batchSize > 0);
        setLearnRate(learnRate_);
    }

    template<Scalar T>
    inline void SGD<T>::recordBegin() {
        to = T(AnyValue);
    }
    
    template<Scalar T>
    inline void SGD<T>::recordEnd() {
        from = T(AnyValue);
    }

    template<Scalar T>
    void SGD<T>::step() const {
        using SegmentType = typename TracerType::SegmentType;
        using DiffScalar = typename TracerType::DiffScalar;
        auto& tracer = TracerType::getInstance();
        tracer.forSegmentInRange(from, to, [this](SegmentType& segment) {
            segment.forNodeInRange(from, to, [this](DiffScalar s) {
                s.setValue(s.getValue() - meanLearnRate * s.getGrad());
            });
        });
    }

    template<Scalar T>
    inline void SGD<T>::zero_grad() const {
        TracerType::getInstance().zero_grad(from, to);
    }

    template<Scalar T>
    void SGD<T>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
        meanLearnRate.swap(obj.meanLearnRate);
        std::swap(batchSize, obj.batchSize);
        from.swap(obj.from);
        to.swap(obj.to);
    }

    template<Scalar T>
    void SGD<T>::setLearnRate(ValueType lr) {
        assert(!lr.isZero() && "[Error]: 0 learn rate does nothing");
        learnRate = lr;
        meanLearnRate = lr / ValueType(batchSize);
    }
}
