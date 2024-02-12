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
        using DiffScalar = typename SegmentType::DiffScalar;
    private:
        std::list<SegmentType> traceList;
    public:
        ~device_obj() = default;
        /* Operations */
        template<class... Args>
        SegmentType& pushSegment(Args&&... args);

        inline void reverse(DiffScalar from, DiffScalar to);
        void reverse_from(DiffScalar from);
        void reverse_to(DiffScalar to);
        void reverse();
        inline void zero_grad(DiffScalar from, DiffScalar to);
        void zero_grad_from(DiffScalar from);
        void zero_grad_to(DiffScalar to);
        void forget(DiffScalar from);

        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func);
        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const;
        /* Getters */
        [[nodiscard]] const std::list<SegmentType>& getTraceList() const noexcept { return traceList; }
        [[nodiscard]] size_t getNumRecord() const noexcept;
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static size_t distance(DiffScalar from, DiffScalar to);
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
    inline void device_obj<DiffTracer<ScalarType>>::reverse(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [](SegmentType& segment) {
            segment.reverse();
        });
    }

    template<class ScalarType>
    void device_obj<DiffTracer<ScalarType>>::reverse_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            segment.reverse();
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void device_obj<DiffTracer<ScalarType>>::reverse_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.reverse();
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void device_obj<DiffTracer<ScalarType>>::reverse() {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite)
            (*ite).reverse();
    }

    template<class ScalarType>
    inline void device_obj<DiffTracer<ScalarType>>::zero_grad(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [](SegmentType& segment) {
            segment.zero_grad();
        });
    }

    template<class ScalarType>
    void device_obj<DiffTracer<ScalarType>>::zero_grad_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;
            segment.zero_grad();
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void device_obj<DiffTracer<ScalarType>>::zero_grad_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.zero_grad();
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void device_obj<DiffTracer<ScalarType>>::forget(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            if (segment.isFound(from)) {
                assert(segment.getLength() == 1 && "[Error]: Forget part of a segment is not supported");
                traceList.erase(ite, end);
                return;
            }
        }
        assert(false && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    template<class Functor>
    void device_obj<DiffTracer<ScalarType>>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false, foundFinalSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            foundFinalSegment |= segment.isFound(to);
            func(std::ref(segment));
            if (foundFinalSegment)
                break;
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    template<class Functor>
    void device_obj<DiffTracer<ScalarType>>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const {
        const_cast<This*>(this)->forSegmentInRange(from, to, func);
    }

    template<class ScalarType>
    size_t device_obj<DiffTracer<ScalarType>>::getNumRecord() const noexcept {
        size_t result = 0;
        for (auto ite = traceList.cbegin(); ite != traceList.cend(); ++ite)
            result += (*ite).getLength();
        return result;
    }

    template<class ScalarType>
    device_obj<DiffTracer<ScalarType>>& device_obj<DiffTracer<ScalarType>>::getInstance() noexcept {
        thread_local static This instance{};
        return instance;
    }

    template<class ScalarType>
    size_t device_obj<DiffTracer<ScalarType>>::distance(DiffScalar from, DiffScalar to) {
        if (&from == &to)
            return 0;

        const auto& tracer = getInstance();
        size_t result = 0;
        tracer.forSegmentInRange(from, to, [from, to, &result](const SegmentType& segment) {
            result += segment.makeFromIndex(from) - segment.makeToIndex(to) + 1;
        });
        return result;
    }
}
