/*
 * Copyright 2023 WeiBo He.
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

#include <list>
#include "TraceSegment.h"

namespace Physica::Core {
    /**
     * \class DiffTracer records a compute graph
     */
    template<class ScalarType>
    class DiffTracer {
        using SegmentType = TraceSegment<ScalarType>;
    public:
        using DiffScalar = typename SegmentType::DiffScalar;
        using DiffRecord = typename SegmentType::DiffRecord;
    private:
        std::list<SegmentType> traceList;
    public:
        ~DiffTracer() = default;
        /* Operations */
        inline void pushOperand(DiffScalar operand);
        template<size_t Size>
        inline void pushOperand(const DiffScalar (&operand)[Size]);

        [[nodiscard]] DiffScalar pushOperation(ScalarType value, ExpressionType source);
        template<size_t Size>
        [[nodiscard]] DiffScalar pushOperation(const SIMD<ScalarType, Size>& simd, ExpressionType source);

        void reverse(DiffScalar from);
        void reverse(DiffScalar from, DiffScalar to);
        inline void zero_grad(DiffScalar from);
        void zero_grad(DiffScalar from, DiffScalar to);
        void forget(DiffScalar from);
        void reserve(size_t size);
        inline void clear();

        [[nodiscard]] inline bool checkLastOpDone() const;
        /* Getters */
        [[nodiscard]] const std::list<SegmentType>& getTraceList() const noexcept { return traceList; }
        [[nodiscard]] size_t getNumRecord() const noexcept;
        [[nodiscard]] ExpressionType getSource(DiffScalar s) const noexcept;
        /* Static members */
        [[nodiscard]] static DiffTracer& getInstance() noexcept;
    private:
        DiffTracer() = default;
        DiffTracer(const DiffTracer&) = default;
        DiffTracer(DiffTracer&&) noexcept = default;
        /* Operators */
        DiffTracer& operator=(const DiffTracer&) = default;
        DiffTracer& operator=(DiffTracer&&) noexcept = default;
        /* Operations */
        void pushSegment(SegmentType segment);
    };

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::pushOperand(DiffScalar operand) {
        assert(!traceList.empty() && !traceList.back().empty() && "[Error]: Pushing operand but operation is unknown");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");
        auto& segment = traceList.back();
        segment.operands.append(operand);
    }

    template<class ScalarType>
    template<size_t Size>
    inline void DiffTracer<ScalarType>::pushOperand(const DiffScalar (&operand)[Size]) {
        assert(!traceList.empty() && !traceList.back().empty() && "[Error]: Pushing operand but operation is unknown");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");

        auto& segment = traceList.back();
        auto& operands = segment.operands;
        const size_t length = operands.getLength();
        const size_t newLength = length + Size;
        operands.reserve(newLength);
        operands.setLength(newLength);

        auto* const p = operands.data() + length;
        for (size_t i = 0; i < Size; ++i)
            p[i] = operand[i];
    }

    template<class ScalarType>
    inline typename DiffTracer<ScalarType>::DiffScalar
    DiffTracer<ScalarType>::pushOperation(ScalarType value, ExpressionType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        if (traceList.empty() || traceList.back().full())
            pushSegment(SegmentType{});
        auto& segment = traceList.back();
        auto& operands = segment.operands;
        segment.records.grow(DiffRecord{operands.getLength(), source});

        auto& values = segment.values;
        const size_t offset = values.getLength();
        values.grow(std::move(value));

        auto& tangents = segment.tangents;
        tangents.grow(0);
        return DiffScalar(values.data() + offset, tangents.data() + offset);
    }

    template<class ScalarType>
    template<size_t Size>
    typename DiffTracer<ScalarType>::DiffScalar
    DiffTracer<ScalarType>::pushOperation(const SIMD<ScalarType, Size>& simd, ExpressionType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        /* Allocate */ {
            if (traceList.empty())
                pushSegment(SegmentType{});
            auto& segment = traceList.back();
            const bool isSpaceEnough = segment.getLength() + Size <= segment.getCapacity();
            if (!isSpaceEnough)
                pushSegment(SegmentType{});
        }

        auto& segment = traceList.back();
        const size_t oldNumRecord = segment.getLength();
        const size_t newNumRecord = oldNumRecord + Size;
        {
            auto& records = segment.records;
            records.setLength(newNumRecord);

            auto& operand = segment.operands;
            const size_t idFirstOperand = operand.getLength();
            const unsigned int numOp = (source == ExpressionType::MulAdd2 || source == ExpressionType::MulAdd4 || source == ExpressionType::MulAdd8) // Until now, these operantions are optimized using SIMD
                                     ? 0
                                     : SegmentType::numOperand(source);
            for (size_t i = 0; i < Size; ++i)
                records[oldNumRecord + i] = DiffRecord{idFirstOperand + i * numOp, source};
        }
        auto& values = segment.values;
        auto* const pValue = values.data_ptr(oldNumRecord);
        simd.store(pValue);
        values.setLength(newNumRecord);

        auto& tangents = segment.tangents;
        SIMD<ScalarType, Size> empty(0);
        auto* const pTangent = tangents.data_ptr(oldNumRecord);
        empty.store(pTangent);
        tangents.setLength(newNumRecord);
        return DiffScalar(pValue, pTangent);
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            segment.reverse(from);
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse(DiffScalar from, DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false, foundFinalSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            foundFinalSegment |= segment.isFound(to);
            segment.reverse(from, to);
            if (foundFinalSegment)
                break;
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::zero_grad(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;
            segment.zero_grad(from);
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::zero_grad(DiffScalar from, DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false, foundFinalSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            foundFinalSegment |= segment.isFound(to);
            segment.zero_grad(from, to);
            if (foundFinalSegment)
                break;
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::forget(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            if (segment.isFound(from)) {
                segment.forget(from);
                traceList.erase(++ite, end);
                return;
            }
        }
        assert(false && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }
    /**
     * Ensure current segment have at least \param size elements of empty space for next continuous store
     */
    template<class ScalarType>
    void DiffTracer<ScalarType>::reserve(size_t size) {
        const auto& segment = traceList.back();
        const bool isSpaceEnough = segment.getLength() + size <= segment.getCapacity();
        if (isSpaceEnough)
            return;
        
        const bool isSmallSize = size < SegmentType::DefaultSize;
        if (isSmallSize)
            pushSegment(SegmentType{});
        else
            pushSegment(SegmentType(size));
    }

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::clear() {
        traceList.clear();
    }

    template<class ScalarType>
    inline bool DiffTracer<ScalarType>::checkLastOpDone() const {
        if (traceList.empty() || traceList.back().empty())
            return true;
        const auto& segment = traceList.back();
        const auto& operands = segment.operands;
        const auto& records = segment.records;
        const auto& record = *records.crbegin();
        return operands.getLength() >= (record.idFirstOperand + SegmentType::numOperand(record.source));
    }

    template<class ScalarType>
    size_t DiffTracer<ScalarType>::getNumRecord() const noexcept {
        size_t result = 0;
        for (auto ite = traceList.cbegin(); ite != traceList.cend(); ++ite)
            result += (*ite).getLength();
        return result;
    }

    template<class ScalarType>
    ExpressionType DiffTracer<ScalarType>::getSource(DiffScalar s) const noexcept {
        for (auto ite = traceList.cbegin(); ite != traceList.cend(); ++ite) {
            const auto& segment = *ite;
            if (!segment.isFound(s))
                continue;
            const size_t index = s.value_ptr() - segment.values().data();
            return segment.records[index].source;
        }
        assert(false && "[Error]: Cannot find the record, maybe it is on another thread?");
        return ExpressionType::Set;
    }

    template<class ScalarType>
    DiffTracer<ScalarType>& DiffTracer<ScalarType>::getInstance() noexcept {
        thread_local static DiffTracer<ScalarType> instance{};
        return instance;
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::pushSegment(SegmentType segment) {
        if (!traceList.empty())
            traceList.back().squeeze();
        traceList.emplace_back(std::move(segment));
    }
}
