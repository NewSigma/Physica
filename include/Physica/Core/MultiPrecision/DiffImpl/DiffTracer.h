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

#include <forward_list>
#include "TraceSegment.h"
#include "Print/ColorfulTracer.h"

namespace Physica::Core {
    /**
     * \class DiffTracer records a compute graph
     */
    template<Scalar T, int Order>
    class DiffTracer {
        static_assert(!T::isDifferentiable, "[Error]: Diff<> pack is not necessary");
        static_assert(Order > 0, "[Error]: 0 order is not differentiable");
        using This = DiffTracer<T, Order>;
    public:
        using SegmentType = TraceSegment<T, Order>;
        using TraceListType = std::forward_list<SegmentType>;
        using DiffScalar = typename SegmentType::DiffScalar;
        using DiffRecord = typename SegmentType::DiffRecord;
    private:
        using ColorfulType = ColorfulTracer<T, Order>;

        TraceListType traceList;
    public:
        ~DiffTracer() = default;
        /* Operators */
        template<Scalar U, unsigned int AnyOrder>
        friend std::ostream& operator<<(std::ostream& os, const DiffTracer<U, AnyOrder>& tracer);
        /* Operations */
        template<class... Args>
        inline void pushOperand(DiffScalar operand, Args... args);
        template<size_t Size>
        inline void pushOperand(const DiffScalar (&operand)[Size]);

        [[nodiscard]] DiffScalar pushOperation(T value, ExprType source);
        template<size_t Size>
        [[nodiscard]] DiffScalar pushOperation(const SIMD<T, Size>& simd, ExprType source);
        [[nodiscard]] DiffScalar copy(DiffScalar s);

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
        void forget(DiffScalar from, DiffScalar to);
        void forget_to(DiffScalar to);
        void reserve(size_t size);
        inline void clear();

        [[nodiscard]] ColorfulType color() const noexcept { return ColorfulType(*this); }

        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func);
        template<class Functor> void forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const;
        [[nodiscard]] inline bool checkLastOpDone() const;
        /* Getters */
        [[nodiscard]] const TraceListType& getTraceList() const noexcept { return traceList; }
        [[nodiscard]] bool isFound(DiffScalar s) const noexcept;
        [[nodiscard]] size_t getNumRecord() const noexcept;
        [[nodiscard]] ExprType getSource(DiffScalar s) const noexcept;
        [[nodiscard]] DiffScalar* getFirstOperand(DiffScalar s);
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
        [[nodiscard]] static DiffRecord makeSetOpNode(const Array<DiffRecord>& records);
    };

    template<Scalar U, unsigned int AnyOrder>
    std::ostream& operator<<(std::ostream& os, const DiffTracer<U, AnyOrder>& tracer) {
        const auto& list = tracer.traceList;
        for (auto ite = list.cbegin(); ite != list.cend(); ++ite) {
            const auto& segment = *ite;
            os << "Trace " << segment.getValues().data() << ":\n";
            os << segment << '\n';
        }
        return os;
    }

    template<Scalar T, int Order>
    template<class... Args>
    inline void DiffTracer<T, Order>::pushOperand(DiffScalar operand, Args... args) {
        assert(!traceList.empty() && !traceList.front().empty() && "[Error]: Pushing operand but operation is unknown");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");
        auto& segment = traceList.front();
        auto& operands = segment.operands;
        const size_t length = operands.getLength();
        const size_t newLength = length + 1 + sizeof...(Args);
        operands.reserve(newLength);
        operands.setLength(newLength);
        pushOperandImpl(operands.data() + length, operand, args...);
    }

    template<Scalar T, int Order>
    template<size_t Size>
    inline void DiffTracer<T, Order>::pushOperand(const DiffScalar (&operand)[Size]) {
        assert(!traceList.empty() && !traceList.front().empty() && "[Error]: Pushing operand but operation is unknown");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");

        auto& segment = traceList.front();
        auto& operands = segment.operands;
        const size_t length = operands.getLength();
        const size_t newLength = length + Size;
        operands.reserve(newLength);
        operands.setLength(newLength);

        auto* const p = operands.data() + length;
        for (size_t i = 0; i < Size; ++i)
            p[i] = operand[i];
    }

    template<Scalar T, int Order>
    inline typename DiffTracer<T, Order>::DiffScalar
    DiffTracer<T, Order>::pushOperation(T value, ExprType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        if (traceList.empty() || traceList.front().full())
            pushSegment();
        auto& segment = traceList.front();
        auto& records = segment.records;
        if (source == ExprType::Set)
            records.grow(makeSetOpNode(records));
        else
            records.grow(DiffRecord{segment.operands.getLength(), source});

        auto& values = segment.values;
        const size_t offset = values.getLength();
        values.grow(value.getValue());

        if constexpr (Order == 1) {
            auto& grads = segment.grads;
            grads.grow(0);
            return DiffScalar(values.data() + offset, grads.data() + offset);
        }
        else {
            auto& tracer = DiffTracer<T, Order - 1>::getInstance();
            auto grad = tracer.pushOperation(T(0), ExprType::Diff);
            segment.grads.grow(grad);
            tracer.pushOperand(grad);
            return DiffScalar(values.data() + offset, std::move(grad));
        }
    }

    template<Scalar T, int Order>
    template<size_t Size>
    typename DiffTracer<T, Order>::DiffScalar
    DiffTracer<T, Order>::pushOperation(const SIMD<T, Size>& simd, ExprType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        /* Allocate */ {
            if (traceList.empty())
                pushSegment();
            auto& segment = traceList.front();
            const bool isSpaceEnough = segment.getLength() + Size <= segment.getCapacity();
            if (!isSpaceEnough)
                pushSegment();
        }

        auto& segment = traceList.front();
        const size_t oldNumRecord = segment.getLength();
        const size_t newNumRecord = oldNumRecord + Size;
        auto& records = segment.records;
        if (source == ExprType::Set) {
            const DiffRecord record = makeSetOpNode(records);
            records.setLength(newNumRecord);
            for (size_t i = 0; i < Size; ++i)
                records[oldNumRecord + i] = record;
        }
        else {
            records.setLength(newNumRecord);
            const size_t idFirstOperand = segment.operands.getLength();
            const unsigned int numOp = (source == ExprType::MulAdd2 || source == ExprType::MulAdd4 || source == ExprType::MulAdd8) // Until now, these operantions are optimized using SIMD
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
        SIMD<T, Size> empty(0);
        auto* const pGrad = grads.data_ptr(oldNumRecord);
        empty.store(pGrad);
        grads.setLength(newNumRecord);
        return DiffScalar(pValue, pGrad);
    }
    /**
     * Usually have seperated nodes continuous on memory to apply SIMD, or communicate between \class DiffTracers.
     */
    template<Scalar T, int Order>
    typename DiffTracer<T, Order>::DiffScalar DiffTracer<T, Order>::copy(DiffScalar s) {
        if constexpr (Order == 1) {
            const auto result = pushOperation(s.getValue(), ExprType::Assign);
            pushOperand(s);
            return result;
        }
        else {
            assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
            assert(s.getSource() == ExprType::Diff && "[Error]: Unexpected source");
            if (traceList.empty() || traceList.front().full())
                pushSegment();

            auto& segment = traceList.front();
            auto& records = segment.records;
            records.grow(makeSetOpNode(records));

            auto& values = segment.values;
            const size_t offset = values.getLength();
            values.grow(s.getValue());

            auto grad = s.getGrad();
            segment.grads.grow(grad);
            return DiffScalar(values.data() + offset, std::move(grad));
        }
    }

    template<Scalar T, int Order>
    template<class... Args>
    typename DiffTracer<T, Order>::SegmentType& DiffTracer<T, Order>::pushSegment(Args&&... args) {
        if (!traceList.empty())
            traceList.front().squeeze();
        return traceList.emplace_front(SegmentType(std::forward<Args>(args)...));
    }

    template<Scalar T, int Order>
    inline void DiffTracer<T, Order>::reverse(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [from, to](SegmentType& segment) {
            segment.reverse(from, to);
        });
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::reverse_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;

            segment.reverse_from(from);
        }
        assert(foundBeginSegment && "[Error]: Cannot find record, maybe it is on another thread?");
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::reverse_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto rend = traceList.end();
        for (auto ite = traceList.begin(); ite != rend; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.reverse_to(to);
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find record, maybe it is on another thread?");
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::reverse() {
        assert(!traceList.empty() && "[Error]: No record found");
        const auto rend = traceList.end();
        for (auto ite = traceList.begin(); ite != rend; ++ite)
            (*ite).reverse();
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::zero_grad(DiffScalar from, DiffScalar to) {
        forSegmentInRange(from, to, [from, to](SegmentType& segment) {
            segment.zero_grad(from, to);
        });
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::zero_grad_from(DiffScalar from) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundBeginSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundBeginSegment |= segment.isFound(from);
            if (!foundBeginSegment)
                continue;
            segment.zero_grad_from(from);
        }
        assert(foundBeginSegment && "[Error]: Cannot find record, maybe it is on another thread?");
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::zero_grad_to(DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        bool foundFinalSegment = false;
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite) {
            auto& segment = *ite;
            foundFinalSegment |= segment.isFound(to);
            segment.zero_grad_to(to);
            if (foundFinalSegment)
                break;
        }
        assert(foundFinalSegment && "[Error]: Cannot find record, maybe it is on another thread?");
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::zero_grad() {
        const auto end = traceList.end();
        for (auto ite = traceList.begin(); ite != end; ++ite)
            (*ite).zero_grad();
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::forget(DiffScalar from, DiffScalar to) {
        assert(!traceList.empty() && "[Error]: No record found");
        auto ite = traceList.begin();
        /* Find from */ {
            const auto end = traceList.end();
            for (; ite != end; ++ite) {
                const bool isFound = (*ite).isFound(from);
                if (isFound)
                    break;
            }
            assert(ite != end && "[Error]: Cannot find record, maybe it is on another thread?");   
        }
        auto ite0 = ite;
        /* Find from */ {
            assert((*ite).getLength() == 1 && "[Error]: Part forget leads to memory fragment");
            ++ite;
            while (ite != traceList.end() && !(*ite).isFound(to))
                ite = traceList.erase_after(ite0);
            assert(ite != traceList.end() && "[Error]: Cannot find record, maybe it is on another thread?");
        }
        auto& s = *ite;
        s.forget_to(to);
        if (s.empty())
            traceList.erase_after(ite0);
    }

    template<Scalar T, int Order>
    void DiffTracer<T, Order>::forget_to(DiffScalar to) {
        while (!traceList.empty() && !traceList.front().isFound(to))
            traceList.pop_front();
        assert(!traceList.empty() && "[Error]: Cannot find record, maybe it is on another thread?");

        auto& front = traceList.front();
        front.forget_to(to);
        if (front.empty())
            traceList.pop_front();
    }
    /**
     * Ensure current segment have at least \param size elements of empty space for next continuous store
     */
    template<Scalar T, int Order>
    void DiffTracer<T, Order>::reserve(size_t size) {
        if (!traceList.empty()) {
            const auto& segment = traceList.front();
            const bool isSpaceEnough = segment.getLength() + size <= segment.getCapacity();
            if (isSpaceEnough)
                return;
        }

        const bool isSmallSize = size < SegmentType::DefaultSize;
        if (isSmallSize)
            pushSegment();
        else
            pushSegment(size);
    }

    template<Scalar T, int Order>
    inline void DiffTracer<T, Order>::clear() {
        traceList.clear();
    }

    template<Scalar T, int Order>
    template<class Functor>
    void DiffTracer<T, Order>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) {
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
        assert(foundBeginSegment && "[Error]: Cannot find record, maybe it is on another thread?");
        assert(foundFinalSegment && "[Error]: Cannot find record, maybe it is on another thread?");
    }

    template<Scalar T, int Order>
    template<class Functor>
    void DiffTracer<T, Order>::forSegmentInRange(DiffScalar from, DiffScalar to, Functor func) const {
        const_cast<This*>(this)->forSegmentInRange(from, to, func);
    }

    template<Scalar T, int Order>
    inline bool DiffTracer<T, Order>::checkLastOpDone() const {
        if (traceList.empty() || traceList.front().empty())
            return true;
        const auto& segment = traceList.front();
        const auto& operands = segment.operands;
        const auto& records = segment.records;
        const auto& record = *records.crbegin();
        return record.source == ExprType::Set
            || operands.getLength() >= (record.idFirstOperand + SegmentType::numOperand(record.source));
    }

    template<Scalar T, int Order>
    bool DiffTracer<T, Order>::isFound(DiffScalar s) const noexcept {
        for (auto& trace : traceList)
            if (trace.isFound(s))
                return true;
        return false;
    }

    template<Scalar T, int Order>
    size_t DiffTracer<T, Order>::getNumRecord() const noexcept {
        size_t result = 0;
        for (auto ite = traceList.cbegin(); ite != traceList.cend(); ++ite)
            result += (*ite).getLength();
        return result;
    }

    template<Scalar T, int Order>
    ExprType DiffTracer<T, Order>::getSource(DiffScalar s) const noexcept {
        for (auto ite = traceList.cbegin(); ite != traceList.cend(); ++ite) {
            const auto& segment = *ite;
            if (!segment.isFound(s))
                continue;
            const size_t index = s.value_ptr() - segment.values.data();
            return segment.records[index].source;
        }
        assert(false && "[Error]: Cannot find the record, maybe it is on another thread?");
        return ExprType::Set;
    }

    template<Scalar T, int Order>
    typename DiffTracer<T, Order>::DiffScalar* DiffTracer<T, Order>::getFirstOperand(DiffScalar s) {
        for (auto ite = traceList.begin(); ite != traceList.end(); ++ite) {
            auto& segment = *ite;
            if (!segment.isFound(s))
                continue;
            const size_t index = s.value_ptr() - segment.values.data();
            return segment.operands.data() + segment.records[index].idFirstOperand;
        }
        assert(false && "[Error]: Cannot find the record, maybe it is on another thread?");
        return {};
    }

    template<Scalar T, int Order>
    DiffTracer<T, Order>& DiffTracer<T, Order>::getInstance() noexcept {
        thread_local static DiffTracer<T, Order> instance{};
        return instance;
    }

    template<Scalar T, int Order>
    size_t DiffTracer<T, Order>::distance(DiffScalar from, DiffScalar to) {
        if (from == to)
            return 0;

        const auto& tracer = getInstance();
        size_t result = 0;
        tracer.forSegmentInRange(from, to, [from, to, &result](const SegmentType& segment) {
            result += segment.makeFromIndex(from) - segment.makeToIndex(to) + 1;
        });
        return result;
    }

    template<Scalar T, int Order>
    template<class... Args>
    inline void DiffTracer<T, Order>::pushOperandImpl(DiffScalar* p, DiffScalar operand, Args... args) {
        *p = operand;
        pushOperandImpl(p + 1, args...);
    }

    template<Scalar T, int Order>
    typename DiffTracer<T, Order>::DiffRecord DiffTracer<T, Order>::makeSetOpNode(const Array<DiffRecord>& records) {
        const size_t length = records.getLength();
        size_t lastSetOpNode = 0;
        if (length != 0) {
            const auto lastRecord = records[length - 1];
            if (lastRecord.source == ExprType::Set)
                lastSetOpNode = lastRecord.idFirstOperand;
            else
                lastSetOpNode = length;
        }
        return DiffRecord{lastSetOpNode, ExprType::Set};
    }
}
