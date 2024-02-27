/*
 * Copyright 2024 WeiBo He.
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

#include "SGD.h"

namespace Physica::Core {
    template<class ScalarType>
    class SGD<device_obj<ScalarType>> {
        static_assert(ScalarType::isDifferentiable, "[Error]: ScalarType must be differentiable");
        using This = SGD<device_obj<ScalarType>>;
        using DeviceScalar = device_obj<ScalarType>;
    public:
        using PlainScalar = typename ScalarType::PlainScalar;
        constexpr static size_t MaxNumThreadPerBlock = 256;
    protected:
        PlainScalar learnRate;
        PlainScalar meanLearnRate;
        unsigned int batchSize;
        DeviceScalar from;
        DeviceScalar to;
    private:
        using TracerType = typename DeviceScalar::TracerType;
        using SegmentType = typename TracerType::SegmentType;
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
        void step();
        void swap(SGD& __restrict obj) noexcept;

        __device__ void stepKernelImpl(SegmentType& segment) const;
        /* Getters */
        [[nodiscard]] PlainScalar getLearnRate() const noexcept { return learnRate; }
        [[nodiscard]] PlainScalar getMeanLearnRate() const noexcept { return meanLearnRate; }
        [[nodiscard]] unsigned int getBatchSize() const noexcept { return batchSize; }
        /* Setters */
        void setLearnRate(PlainScalar lr);
    };
}

#include "SGDImpl.cuh"
