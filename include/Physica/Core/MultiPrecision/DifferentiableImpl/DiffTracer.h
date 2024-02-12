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

#include <list>
#include "TraceSegment.h"

namespace Physica::Core {
    /**
     * \class DiffTracer records a compute graph
     */
    template<class ScalarType>
    class DiffTracer {
        using This = DiffTracer<ScalarType>;
    public:
        using SegmentType = TraceSegment<ScalarType>;
        using DiffScalar = typename SegmentType::DiffScalar;
        using DiffRecord = typename SegmentType::DiffRecord;
    private:
        std::list<SegmentType> traceList;
    public:
        ~DiffTracer() = default;
        /* Operations */
        template<class... Args>
        inline void pushOperand(DiffScalar operand, Args... args);
        template<size_t Size>
        inline void pushOperand(const DiffScalar (&operand)[Size]);

        [[nodiscard]] DiffScalar pushOperation(ScalarType value, ExpressionType source);
        template<size_t Size>
        [[nodiscard]] DiffScalar pushOperation(const SIMD<ScalarType, Size>& simd, ExpressionType source);

        template<class... Args>
        SegmentType& pushSegment(Args&&... args);

        inline void reverse(DiffScalar from, DiffScalar to);
        void reverse_from(DiffScalar from);
        void reverse_to(DiffScalar to);
        void reverse();
        inline void zero_grad(DiffScalar from, DiffScalar to);
        void zero_grad_from(DiffScalar from);
        void zero_grad_to(DiffScalar to);
        void zero_grad();
        void forget(DiffScalar from);
        void reserve(size_t size);
        inline void clear();

        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func);
        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const;
        [[nodiscard]] inline bool checkLastOpDone() const;
        /* Getters */
        [[nodiscard]] const std::list<SegmentType>& getTraceList() const noexcept { return traceList; }
        [[nodiscard]] size_t getNumRecord() const noexcept;
        [[nodiscard]] ExpressionType getSource(DiffScalar s) const noexcept;
        /* Static members */
        [[nodiscard]] static DiffTracer& getInstance() noexcept;
        [[nodiscard]] static size_t distance(DiffScalar from, DiffScalar to);
    private:
        DiffTracer() = default;
        DiffTracer(const DiffTracer&) = default;
        DiffTracer(DiffTracer&&) noexcept = default;
        /* Operators */
        DiffTracer& operator=(const DiffTracer&) = default;
        DiffTracer& operator=(DiffTracer&&) noexcept = default;
        /* Operations */
        template<class... Args>
        inline void pushOperandImpl(DiffScalar* p, DiffScalar operand, Args... args);
        void pushOperandImpl([[maybe_unused]] DiffScalar* p) {}
        /* Static members */
        [[nodiscard]] static DiffRecord makeSetOpNode(const Utils::Array<DiffRecord>& records);
    };

    template<class ScalarType>
    template<class... Args>
    inline void DiffTracer<ScalarType>::pushOperand(DiffScalar operand, Args... args) {
        assert(!traceList.empty() && !traceList.back().empty() && "[Error]: Pushing operand but operation is unknown");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");
        auto& segment = traceList.back();
        auto& operands = segment.operands;
        const size_t length = operands.getLength();
        const size_t newLength = length + 1 + sizeof...(Args);
        operands.reserve(newLength);
        operands.setLength(newLength);
        pushOperandImpl(operands.data() + length, operand, args...);
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
            pushSegment();
        auto& segment = traceList.back();
        auto& records = segment.records;
        if (source == ExpressionType::Set)
            records.grow(makeSetOpNode(records));
        else
            records.grow(DiffRecord{segment.operands.getLength(), source});

        auto& values = segment.values;
        const size_t offset = values.getLength();
        values.grow(std::move(value));

        auto& grads = segment.grads;
        grads.grow(0);
        return DiffScalar(values.data() + offset, grads.data() + offset);
    }

    template<class ScalarType>
    template<size_t Size>
    typename DiffTracer<ScalarType>::DiffScalar
    DiffTracer<ScalarType>::pushOperation(const SIMD<ScalarType, Size>& simd, ExpressionType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        /* Allocate */ {
            if (traceList.empty())
                pushSegment();
            auto& segment = traceList.back();
            const bool isSpaceEnough = segment.getLength() + Size <= segment.getCapacity();
            if (!isSpaceEnough)
                pushSegment();
        }

        auto& segment = traceList.back();
        const size_t oldNumRecord = segment.getLength();
        const size_t newNumRecord = oldNumRecord + Size;
        auto& records = segment.records;
        if (source == ExpressionType::Set) {
            const DiffRecord record = makeSetOpNode(records);
            records.setLength(newNumRecord);
            for (size_t i = 0; i < Size; ++i)
                records[oldNumRecord + i] = record;
        }
        else {
            records.setLength(newNumRecord);
            const size_t idFirstOperand = segment.operands.getLength();
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

        auto& grads = segment.grads;
        SIMD<ScalarType, Size> empty(0);
        auto* const pGrad = grads.data_ptr(oldNumRecord);
        empty.store(pGrad);
        grads.setLength(newNumRecord);
        return DiffScalar(pValue, pGrad);
    }

    template<class ScalarType>
    template<class... Args>
    typename DiffTracer<ScalarType>::SegmentType& DiffTracer<ScalarType>::pushSegment(Args&&... args) {
        if (!traceList.empty())
            traceList.back().squeeze();
        return traceList.emplace_back(SegmentType(std::forward<Args>(args)...));
    }

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::reverse(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [from, to](SegmentType& segment) {
            segment.reverse(from, to);
        });
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            segment.reverse_from(from);
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.reverse_to(to);
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse() {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite)
            (*ite).reverse();
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::zero_grad(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [from, to](SegmentType& segment) {
            segment.zero_grad(from, to);
        });
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::zero_grad_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;
            segment.zero_grad_from(from);
        }
        assert(foundBeginSegment && "[Error]: Cannot find begin record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::zero_grad_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.zero_grad_to(to);
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find final record, maybe it is on another thread?");
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::zero_grad() {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto rend = traceList.rend();
        for (auto ite = traceList.rbegin(); ite != rend; ++ite)
            (*ite).zero_grad();
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
            pushSegment();
        else
            pushSegment(size);
    }

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::clear() {
        traceList.clear();
    }

    template<class ScalarType>
    template<class Functor>
    void DiffTracer<ScalarType>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) {
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
    void DiffTracer<ScalarType>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const {
        const_cast<This*>(this)->forSegmentInRange(from, to, func);
    }

    template<class ScalarType>
    inline bool DiffTracer<ScalarType>::checkLastOpDone() const {
        if (traceList.empty() || traceList.back().empty())
            return true;
        const auto& segment = traceList.back();
        const auto& operands = segment.operands;
        const auto& records = segment.records;
        const auto& record = *records.crbegin();
        return record.source == ExpressionType::Set
            || operands.getLength() >= (record.idFirstOperand + SegmentType::numOperand(record.source));
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
            const size_t index = s.value_ptr() - segment.values.data();
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
    size_t DiffTracer<ScalarType>::distance(DiffScalar from, DiffScalar to) {
        if (&from == &to)
            return 0;

        const auto& tracer = getInstance();
        size_t result = 0;
        tracer.forSegmentInRange(from, to, [from, to, &result](const SegmentType& segment) {
            result += segment.makeFromIndex(from) - segment.makeToIndex(to) + 1;
        });
        return result;
    }

    template<class ScalarType>
    template<class... Args>
    inline void DiffTracer<ScalarType>::pushOperandImpl(DiffScalar* p, DiffScalar operand, Args... args) {
        *p = operand;
        pushOperandImpl(p + 1, args...);
    }

    template<class ScalarType>
    typename DiffTracer<ScalarType>::DiffRecord DiffTracer<ScalarType>::makeSetOpNode(const Utils::Array<DiffRecord>& records) {
        const size_t length = records.getLength();
        size_t lastSetOpNode = 0;
        if (length != 0) {
            const auto lastRecord = records[length - 1];
            if (lastRecord.source == ExpressionType::Set)
                lastSetOpNode = lastRecord.idFirstOperand;
            else
                lastSetOpNode = length;
        }
        return DiffRecord{lastSetOpNode, ExpressionType::Set};
    }
}
