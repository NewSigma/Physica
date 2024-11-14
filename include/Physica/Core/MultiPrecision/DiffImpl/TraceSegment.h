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

#include "Physica/Core/MultiPrecision/ExprType.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Print/ColorfulTraceSegment.h"

namespace Physica::Core {
    template<class ScalarType, int Order> class DiffTracer;
    /**
     * \class TraceSegment provides a segment of continuous memory to record compute graph
     */
    template<class ScalarType, int Order>
    class TraceSegment {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Diff<> pack is not necessary");
        static_assert(Order > 0, "[Error]: 0 order is not differentiable");
        using This = TraceSegment<ScalarType, Order>;
        using GradType = typename std::conditional<Order == 1, ScalarType, Diff<ScalarType, DiffMode::Reverse, Order - 1>>::type;
    public:
        constexpr static size_t DefaultSize = 4096; // I guess it is not a bad choice
        using DiffScalar = Diff<ScalarType, DiffMode::Reverse, Order>;
        using ValueVector = Vector<ScalarType>;
        using GradVector = Vector<GradType>;
        struct DiffRecord {
            using device_obj_type = DiffRecord;

            size_t idFirstOperand;
            ExprType source;
        };
    private:
        using ColorfulType = ColorfulTraceSegment<ScalarType, Order>;

        Array<DiffRecord> records;
        Array<DiffScalar> operands;
        ValueVector values;
        GradVector grads;
    public:
        explicit TraceSegment(size_t size = DefaultSize);
        TraceSegment(ValueVector values_);
        TraceSegment(const TraceSegment&) = default;
        TraceSegment(TraceSegment&&) noexcept = default;
        ~TraceSegment() = default;
        /* Operators */
        TraceSegment& operator=(TraceSegment obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] inline DiffScalar operator[](size_t index);
        [[nodiscard]] inline const DiffScalar operator[](size_t index) const;
        template<Scalar T, unsigned int AnyOrder>
        friend std::ostream& operator<<(std::ostream& os, const TraceSegment<T, AnyOrder>& segment);
        /* Operations */
        inline void reverse(DiffScalar from, DiffScalar to);
        void reverse_from(DiffScalar from) { reverse(makeFromIndex(from), 0); }
        void reverse_to(DiffScalar to) { reverse(getLength() - 1, makeToIndex(to)); }
        inline void reverse();
        void zero_grad(DiffScalar from, DiffScalar to);
        inline void zero_grad_from(DiffScalar from);
        inline void zero_grad_to(DiffScalar to);
        inline void zero_grad();
        void forget_to(DiffScalar to);
        void squeeze();

        [[nodiscard]] ColorfulType color() const noexcept { return ColorfulType(*this); }

        template<class Functor> inline void forNodeInRange(DiffScalar from, DiffScalar to, Functor func);
        template<class Functor> inline void forNodeInRange(DiffScalar from, DiffScalar to, Functor func) const;
        void swap(TraceSegment& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return records.getLength(); }
        [[nodiscard]] size_t getCapacity() const noexcept { return records.getCapacity(); }
        [[nodiscard]] bool empty() const noexcept { return records.empty(); }
        [[nodiscard]] bool full() const noexcept { return records.full(); }
        [[nodiscard]] bool isFound(DiffScalar s) const noexcept { return find(s) < getLength(); }
        [[nodiscard]] size_t find(DiffScalar s) const noexcept;
        [[nodiscard]] const Array<DiffRecord>& getRecords() const noexcept { return records; }
        [[nodiscard]] const Array<DiffScalar>& getOperands() const noexcept { return operands; }
        [[nodiscard]] ValueVector& getValues() noexcept { return values; }
        [[nodiscard]] const ValueVector& getValues() const noexcept { return values; }
        [[nodiscard]] GradVector& getGrads() noexcept { return grads; }
        [[nodiscard]] const GradVector& getGrads() const noexcept { return grads; }
        /* Static members */
        [[nodiscard]] constexpr static unsigned int numOperand(ExprType type);
    private:
        /* Operations */
        void reverse(size_t fromIndex, size_t toIndex);
        template<size_t Size>
        void reverseMulAdd(DiffScalar firstOpX, DiffScalar firstOpY, DiffScalar firstOpZ, size_t traceId);
        [[nodiscard]] size_t makeFromIndex(DiffScalar from) const noexcept;
        [[nodiscard]] size_t makeToIndex(DiffScalar to) const noexcept;
        /* Static members */
        inline static void updateGrad(GradType& target, GradType deltaGrad);
        /* Friends */
        friend class DiffTracer<ScalarType, Order>;
        friend class device_obj<This>;
    };

    template<class ScalarType, int Order>
    TraceSegment<ScalarType, Order>::TraceSegment(size_t size) {
        assert(size > 0 && "[Error]: Bad size");
        records.reserve(size);
        operands.reserve(3 * size); //MulAdd operation is 3-operand
        values.reserve(size);
        grads.reserve(size);
    }

    template<class ScalarType, int Order>
    TraceSegment<ScalarType, Order>::TraceSegment(ValueVector values_) : values(std::move(values_)) {
        const size_t length = values.getLength();
        assert(length > 0 && "[Error]: Bad size");
        records.resize(length, DiffRecord{0, ExprType::Set});
        grads.resize(length, ScalarType(0));
    }

    template<class ScalarType, int Order>
    inline typename TraceSegment<ScalarType, Order>::DiffScalar TraceSegment<ScalarType, Order>::operator[](size_t index) {
        if constexpr (Order == 1)
            return DiffScalar(values.data_ptr(index), grads.data_ptr(index));
        else
            return DiffScalar(values.data_ptr(index), grads[index]);
    }

    template<class ScalarType, int Order>
    inline const typename TraceSegment<ScalarType, Order>::DiffScalar TraceSegment<ScalarType, Order>::operator[](size_t index) const {
        return const_cast<This&>(*this).operator[](index);
    }

    template<Scalar T, unsigned int AnyOrder>
    std::ostream& operator<<(std::ostream& os, const TraceSegment<T, AnyOrder>& segment) {
        const size_t length = segment.getLength();
        for (size_t i = 0; i < length; ++i) {
            const auto source = segment.records[i].source;
            const auto width = os.width();
            os << std::setw(4) << i << std::setw(width) << ": ";
            os << std::setw(10) << exprTypeToStr(source) << std::setw(width) << ' ';
            os << segment.values[i] << ' ' << segment.grads[i] << ' ';
            os << "Op: ";
            const size_t idFirstOperand = segment.records[i].idFirstOperand;
            const size_t num = TraceSegment<T, AnyOrder>::numOperand(source);
            for (size_t j = 0; j < num; ++j) {
                const auto& op = segment.operands[idFirstOperand + j];
                const size_t index = segment.find(op);
                const bool found = index < length;
                if (found)
                    os << index;
                else
                    os << op.value_ptr();
                os << ' ';
            }
            os << '\n';
        }
        return os;
    }

    template<class ScalarType, int Order>
    inline void TraceSegment<ScalarType, Order>::reverse(DiffScalar from, DiffScalar to) {
        assert(!empty());
        reverse(makeFromIndex(from), makeToIndex(to));
    }

    template<class ScalarType, int Order>
    inline void TraceSegment<ScalarType, Order>::reverse() {
        assert(!empty());
        reverse(getLength() - 1, 0);
    }

    template<class ScalarType, int Order>
    void TraceSegment<ScalarType, Order>::zero_grad(DiffScalar from, DiffScalar to) {
        const size_t fromIndex = makeFromIndex(from);
        const size_t toIndex = makeToIndex(to);
        assert(toIndex <= fromIndex && "[Error]: Invalid range");
        auto segment = grads.segment(toIndex, fromIndex + 1);
        if constexpr (Order == 1)
            segment = ScalarType(0);
        else {
            for (size_t i = 0; i < segment.getLength(); ++i)
                segment[i].setValue(0);
        }
    }

    template<class ScalarType, int Order>
    inline void TraceSegment<ScalarType, Order>::zero_grad_from(DiffScalar from) {
        auto segment = grads.segment(0, makeFromIndex(from) + 1);
        if constexpr (Order == 1)
            segment = ScalarType(0);
        else {
            for (size_t i = 0; i < segment.getLength(); ++i)
                segment[i].setValue(0);
        }
    }

    template<class ScalarType, int Order>
    inline void TraceSegment<ScalarType, Order>::zero_grad_to(DiffScalar to) {
        auto segment = grads.segment(makeToIndex(to), getLength());
        if constexpr (Order == 1)
            segment = ScalarType(0);
        else {
            for (size_t i = 0; i < segment.getLength(); ++i)
                segment[i].setValue(0);
        }
    }

    template<class ScalarType, int Order>
    inline void TraceSegment<ScalarType, Order>::zero_grad() {
        if constexpr (Order == 1)
            grads = ScalarType(0);
        else {
            for (auto& grad : grads)
                grad.setValue(0);
        }
    }

    template<class ScalarType, int Order>
    void TraceSegment<ScalarType, Order>::forget_to(DiffScalar to) {
        assert(isFound(to) && "[Error]: forgeting a non existent, this may be a bug");
        const size_t index = find(to);
        const auto& record = records[index];
        if (record.source == ExprType::Set) {
            size_t idLastOperand = 0;
            if (index != 0) {
                const auto& lastRecord = records[index - 1];
                idLastOperand = lastRecord.idFirstOperand + numOperand(lastRecord.source);
            }
            operands.resize(idLastOperand);
        }
        else
            operands.resize(record.idFirstOperand);
        records.resize(index);
        values.resize(index);
        grads.resize(index);
    }
    /**
     * FIXME: squeeze values and grads may invalidate pointers, so the unused memory cannot be freed
     */
    template<class ScalarType, int Order>
    void TraceSegment<ScalarType, Order>::squeeze() {
        records.squeeze();
        operands.squeeze();
    }

    template<class ScalarType, int Order>
    template<class Functor>
    inline void TraceSegment<ScalarType, Order>::forNodeInRange(DiffScalar from, DiffScalar to, Functor func) {
        const size_t to_i = makeToIndex(to);
        const size_t from_i = makeFromIndex(from) + 1;
        for (size_t i = to_i; i < from_i; ++i)
            func(this->operator[](i));
    }

    template<class ScalarType, int Order>
    template<class Functor>
    inline void TraceSegment<ScalarType, Order>::forNodeInRange(DiffScalar from, DiffScalar to, Functor func) const {
        const_cast<This*>(this)->forNodeInRange(from, to, func);
    }

    template<class ScalarType, int Order>
    void TraceSegment<ScalarType, Order>::swap(TraceSegment& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        records.swap(obj.records);
        operands.swap(obj.operands);
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<class ScalarType, int Order>
    size_t TraceSegment<ScalarType, Order>::find(DiffScalar s) const noexcept {
        assert(values.getLength() == grads.getLength() && "[Error]: Invalid param");
        const auto* pValue = values.data();
        const auto* pValue1 = s.value_ptr();
        const size_t length = values.getLength();
        if (pValue1 < pValue)
            return length;
        const size_t index = pValue1 - pValue;
        [[maybe_unused]] const bool isValueFound = index < length;
        assert((!isValueFound || (s == this->operator[](index))) && "[Error]: Bad DiffScalar");
        return index;
    }

    template<class ScalarType, int Order>
    constexpr unsigned int TraceSegment<ScalarType, Order>::numOperand(ExprType type) {
        switch (type) {
            case ExprType::Set: return 0;
            case ExprType::Assign: return 1;
            case ExprType::Diff: return 1;
            case ExprType::Minus: return 1;
            case ExprType::Add: return 2;
            case ExprType::Sub: return 2;
            case ExprType::Mul: return 2;
            case ExprType::Div: return 2;
            case ExprType::Sum: return 2;
            case ExprType::MulAdd2: [[fallthrough]];
            case ExprType::MulAdd4: [[fallthrough]];
            case ExprType::MulAdd8: return 3;
            case ExprType::More: return 2;
            case ExprType::MoreEq: return 2;
            case ExprType::Reciprocal: return 1;
            case ExprType::Sqrt: return 1;
            case ExprType::Cbrt: return 1;
            case ExprType::Abs: return 1;
            case ExprType::Relu: return 1;
            case ExprType::Square: return 1;
            case ExprType::Ln: return 1;
            case ExprType::Exp: return 1;
            case ExprType::Pow: return 2;
            case ExprType::Sin: return 1;
            case ExprType::Cos: return 1;
            case ExprType::ArcCos: return 1;
            case ExprType::Tanh: return 1;
            default: [[unlikely]]
                throw std::invalid_argument("[Error]: Unrecognized type");
        }
    }

    template<class ScalarType, int Order>
    void TraceSegment<ScalarType, Order>::reverse(size_t fromIndex, size_t toIndex) {
        assert(toIndex <= fromIndex && "[Error]: Invalid index");
        for (size_t i = fromIndex; toIndex <= i && i <= fromIndex; --i) {
            const DiffRecord record = records[i];
            const size_t idFirstOperand = record.idFirstOperand;
            const ExprType source = record.source;
            if (source == ExprType::Set) {
                i = idFirstOperand;
                continue;
            }

            DiffScalar currentNode = (*this)[i];
            DiffScalar operandX = operands[idFirstOperand];
            GradType& gradX = operandX.getGrad();
            if (source == ExprType::Diff) {
                auto& grad = currentNode.getGrad();
                if constexpr (Order == 1)
                    grad += gradX;
                else
                    grad.setValue(ScalarType(grad) + ScalarType(gradX));
                continue;
            }

            GradType grad;
            if constexpr (Order == 1)
                grad = currentNode.getGrad();
            else
                grad = *currentNode.getGrad().getFirstOperand();

            if (grad.isZero())
                continue;
            /* Unitary Operations */ {
                switch (source) {
                    case ExprType::Assign:
                    case ExprType::Minus:
                        updateGrad(gradX, grad * GradType(source == ExprType::Assign ? 1.0 : -1.0));
                        continue;
                    case ExprType::Reciprocal:
                        updateGrad(gradX, -grad * square(GradType(currentNode)));
                        continue;
                    case ExprType::Sqrt:
                        updateGrad(gradX, grad / GradType(currentNode) * GradType(0.5));
                        continue;
                    case ExprType::Cbrt:
                        updateGrad(gradX, grad / (square(GradType(currentNode)) * GradType(3)));
                        continue;
                    case ExprType::Abs:
                        updateGrad(gradX, operandX.getValue().isPositive() ? grad : -grad);
                        continue;
                    case ExprType::Relu:
                        updateGrad(gradX, operandX.getValue().isPositive() ? grad : GradType(0));
                        continue;
                    case ExprType::Square:
                        updateGrad(gradX, grad * GradType(operandX) * GradType(2));
                        continue;
                    case ExprType::Ln:
                        updateGrad(gradX, grad / GradType(operandX));
                        continue;
                    case ExprType::Exp:
                        updateGrad(gradX, grad * GradType(currentNode));
                        continue;
                    case ExprType::Sin:
                        updateGrad(gradX, grad * cos(GradType(operandX)));
                        continue;
                    case ExprType::Cos:
                        updateGrad(gradX, -grad * sin(GradType(operandX)));
                        continue;
                    case ExprType::ArcCos:
                        updateGrad(gradX, -grad / sqrt(GradType(1) - square(GradType(operandX))));
                        continue;
                    case ExprType::Tanh:
                        updateGrad(gradX, grad - grad * square(GradType(currentNode)));
                        continue;
                    default:;
                }
            }
            /* Binary Operations */ {
                DiffScalar operandY = operands[idFirstOperand + 1];
                GradType& gradY = operandY.getGrad();
                switch (source) {
                    case ExprType::Add:
                    case ExprType::Sub:
                        updateGrad(gradX, grad);
                        updateGrad(gradY, grad * GradType(source == ExprType::Add ? 1.0 : -1.0));
                        continue;
                    case ExprType::Mul:
                        updateGrad(gradX, grad * GradType(operandY));
                        updateGrad(gradY, grad * GradType(operandX));
                        continue;
                    case ExprType::Div: {
                        const GradType dx = grad * reciprocal(GradType(operandY));
                        updateGrad(gradX, dx);
                        updateGrad(gradY, -dx * GradType(currentNode));
                        continue;
                    }
                    case ExprType::MulAdd2:
                        if constexpr (Order == 1) {
                            if constexpr (ScalarType::Option == Double) {
                                DiffScalar operandZ = operands[idFirstOperand + 2];
                                i -= 1;
                                [[maybe_unused]] const bool isInRange = toIndex <= i && i <= fromIndex;
                                assert(isInRange && "[Error]: Unexpected id");
                                reverseMulAdd<2>(operandX, operandY, operandZ, i);
                            }
                            else {
                                assert(false && "[Error]: MulAdd2 apply to double type only, you should not arrive here");
                            }
                        }
                        else {
                            assert(false && "[Error]: MulAdd2 for high order autodiff not implemented");
                        }
                        continue;
                    case ExprType::MulAdd4: {
                        if constexpr (Order == 1) {
                            DiffScalar operandZ = operands[idFirstOperand + 2];
                            i -= 3;
                            [[maybe_unused]] const bool isInRange = toIndex <= i && i <= fromIndex;
                            assert(isInRange && "[Error]: Unexpected id");
                            reverseMulAdd<4>(operandX, operandY, operandZ, i);
                        }
                        else {
                            assert(false && "[Error]: MulAdd4 for high order autodiff not implemented");
                        }
                        continue;
                    }
                    case ExprType::MulAdd8: {
                        if constexpr (Order == 1) {
                            DiffScalar operandZ = operands[idFirstOperand + 2];
                            i -= 7;
                            [[maybe_unused]] const bool isInRange = toIndex <= i && i <= fromIndex;
                            assert(isInRange && "[Error]: Unexpected id");
                            reverseMulAdd<8>(operandX, operandY, operandZ, i);
                        }
                        else {
                            assert(false && "[Error]: MulAdd8 for high order autodiff not implemented");
                        }
                        continue;
                    }
                    default:;
                }
            }
            assert(false && "[Error]: Undefined operator for back propagation");
        }
    }

    template<class ScalarType, int Order>
    template<size_t Size>
    void TraceSegment<ScalarType, Order>::reverseMulAdd(DiffScalar firstOpX, DiffScalar firstOpY, DiffScalar firstOpZ, size_t traceId) {
        using PacketType = SIMD<ScalarType, Size>;
        const PacketType grad = grads.template packet<PacketType>(traceId);
        PacketType gradZ{};
        gradZ.load(firstOpZ.grad_ptr());
        gradZ += grad;
        gradZ.store(firstOpZ.grad_ptr());


        PacketType valueX{}, valueY{}, gradX{}, gradY{};
        valueX.load(firstOpX.value_ptr());
        valueY.load(firstOpY.value_ptr());
        gradX.load(firstOpX.grad_ptr());
        gradY.load(firstOpY.grad_ptr());
        gradX += grad * valueY;
        gradY += grad * valueX;
        valueX.store(firstOpX.value_ptr());
        valueY.store(firstOpY.value_ptr());
        gradX.store(firstOpX.grad_ptr());
        gradY.store(firstOpY.grad_ptr());
    }

    template<class ScalarType, int Order>
    size_t TraceSegment<ScalarType, Order>::makeFromIndex(DiffScalar from) const noexcept {
        assert(!empty());
        size_t fromIndex = find(from);
        const bool isFromNotFound = fromIndex >= getLength();
        if (isFromNotFound)
            fromIndex = getLength() - 1;
        return fromIndex;
    }

    template<class ScalarType, int Order>
    size_t TraceSegment<ScalarType, Order>::makeToIndex(DiffScalar to) const noexcept {
        assert(!empty());
        size_t toIndex = find(to);
        const bool isToNotFound = toIndex >= getLength();
        if (isToNotFound)
            toIndex = 0;
        return toIndex;
    }

    template<class ScalarType, int Order>
    inline void TraceSegment<ScalarType, Order>::updateGrad(GradType& target, GradType deltaGrad) {
        if constexpr (Order == 1)
            target += deltaGrad;
        else {
            assert(target.getSource() == ExprType::Diff);
            const GradType temp = target + deltaGrad;
            *target.getFirstOperand() = temp;
            target.setValue(ScalarType(temp));
        }
    }
}
