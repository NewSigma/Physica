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

#include "DiffTracer.h"
#include "TraceSegment.cuh"

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<DiffTracer<ScalarType>> {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        using host_obj = DiffTracer<ScalarType>;
        using This = device_obj<host_obj>;
    public:
        using SegmentType = device_obj<typename host_obj::SegmentType>;
    private:
        std::list<SegmentType> traceList;
    public:
        ~device_obj() = default;
        /* Operations */
        template<class... Args>
        SegmentType& pushSegment(Args&&... args);
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
    private:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
    };

    template<class ScalarType>
    template<class... Args>
    typename device_obj<DiffTracer<ScalarType>>::SegmentType& device_obj<DiffTracer<ScalarType>>::pushSegment(Args&&... args) {
        return traceList.emplace_back(SegmentType(std::forward<Args>(args)...));
    }

    template<class ScalarType>
    device_obj<DiffTracer<ScalarType>>& device_obj<DiffTracer<ScalarType>>::getInstance() noexcept {
        thread_local static This instance{};
        return instance;
    }
}
