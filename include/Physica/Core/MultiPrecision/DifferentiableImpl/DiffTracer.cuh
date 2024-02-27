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
    template<class ScalarType, unsigned int Order>
    class device_obj<DiffTracer<ScalarType, Order>> {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        using host_obj = DiffTracer<ScalarType, Order>;
        using This = device_obj<host_obj>;
    public:
        using SegmentType = device_obj<typename host_obj::SegmentType>;
        using TraceListType = std::forward_list<SegmentType>;
        using DiffScalar = typename SegmentType::DiffScalar;
    private:
        using VectorType = typename SegmentType::VectorType;

        TraceListType traceList; //forward_list is FILO
        std::list<SegmentType> reserveList; //list is FIFO
    public:
        ~device_obj() = default;
        /* Operations */
        SegmentType& pushSegment(size_t size, ExpressionType type);
        SegmentType& pushSegment(ScalarType value);
        SegmentType& pushSegment(const VectorType& value);

        inline void reverse(DiffScalar from, DiffScalar to);
        void reverse_from(DiffScalar from);
        void reverse_to(DiffScalar to);
        void reverse();
        inline void zero_grad(DiffScalar from, DiffScalar to);
        void zero_grad_from(DiffScalar from);
        void zero_grad_to(DiffScalar to);
        void forget(DiffScalar to);

        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func);
        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const;
        /* Getters */
        [[nodiscard]] const TraceListType& getTraceList() const noexcept { return traceList; }
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

    template<class ScalarType, unsigned int Order>
    typename device_obj<DiffTracer<ScalarType, Order>>::SegmentType& device_obj<DiffTracer<ScalarType, Order>>::pushSegment(size_t size, ExpressionType type) {
        const auto end = reserveList.end();
        for (auto ite = reserveList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            const bool isSizeMatch = segment.getLength() == size;
            const bool isNumOperandMatch = segment.getNumOperands() == size * SegmentType::numOperand(type);
            if (isSizeMatch && isNumOperandMatch) {
                traceList.push_front(std::move(segment));
                reserveList.erase(ite);
                return traceList.front();
            }
        }

        try {
            return traceList.emplace_front(size, type);
        }
        catch (CudaException& e) {
            if (e.getCode() != cudaErrorMemoryAllocation)
                throw e;
        }
        reserveList.clear();
        return traceList.emplace_front(size, type);
    }

    template<class ScalarType, unsigned int Order>
    typename device_obj<DiffTracer<ScalarType, Order>>::SegmentType& device_obj<DiffTracer<ScalarType, Order>>::pushSegment(ScalarType value) {
        auto& result = pushSegment(1, ExpressionType::Set);
        result.init(value);
        return result;
    }

    template<class ScalarType, unsigned int Order>
    typename device_obj<DiffTracer<ScalarType, Order>>::SegmentType& device_obj<DiffTracer<ScalarType, Order>>::pushSegment(const VectorType& value) {
        auto& result = pushSegment(value.getLength(), ExpressionType::Set);
        result.init(value);
        return result;
    }

    template<class ScalarType, unsigned int Order>
    inline void device_obj<DiffTracer<ScalarType, Order>>::reverse(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [](SegmentType& segment) {
            segment.reverse();
        });
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<DiffTracer<ScalarType, Order>>::reverse_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            segment.reverse();
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<DiffTracer<ScalarType, Order>>::reverse_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.reverse();
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<DiffTracer<ScalarType, Order>>::reverse() {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite)
            (*ite).reverse();
    }

    template<class ScalarType, unsigned int Order>
    inline void device_obj<DiffTracer<ScalarType, Order>>::zero_grad(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [](SegmentType& segment) {
            segment.zero_grad();
        });
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<DiffTracer<ScalarType, Order>>::zero_grad_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;
            segment.zero_grad();
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<DiffTracer<ScalarType, Order>>::zero_grad_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.zero_grad();
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<DiffTracer<ScalarType, Order>>::forget(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ite = traceList.begin()) {
            auto& segment = *ite;
            const bool isFound = segment.isFound(to);
            assert((!isFound || (segment.getLength() == 1)) && "[Error]: Forget part of a segment is not supported");
            reserveList.push_back(std::move(segment));
            traceList.pop_front();
            if (isFound)
                return;
        }
        assert(false && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType, unsigned int Order>
    template<class Functor>
    void device_obj<DiffTracer<ScalarType, Order>>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false, foundFinalSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
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

    template<class ScalarType, unsigned int Order>
    template<class Functor>
    void device_obj<DiffTracer<ScalarType, Order>>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const {
        const_cast<This*>(this)->forSegmentInRange(from, to, func);
    }

    template<class ScalarType, unsigned int Order>
    size_t device_obj<DiffTracer<ScalarType, Order>>::getNumRecord() const noexcept {
        size_t result = 0;
        for (auto ite = traceList.cbegin(); ite != traceList.cend(); ++ite)
            result += (*ite).getLength();
        return result;
    }

    template<class ScalarType, unsigned int Order>
    device_obj<DiffTracer<ScalarType, Order>>& device_obj<DiffTracer<ScalarType, Order>>::getInstance() noexcept {
        thread_local static This instance{};
        return instance;
    }

    template<class ScalarType, unsigned int Order>
    size_t device_obj<DiffTracer<ScalarType, Order>>::distance(DiffScalar from, DiffScalar to) {
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
