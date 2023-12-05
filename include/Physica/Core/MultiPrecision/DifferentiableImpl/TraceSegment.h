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

#include "Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    /**
     * \class TraceSegment provides a segment of continuous memory to record compute graph
     */
    template<class ScalarType>
    class TraceSegment {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
    public:
        constexpr static size_t DefaultSize = 4096; // I guess it is not a bad choice
        using DiffScalar = Differentiable<ScalarType, DiffMode::Reverse>;
        using VectorType = Vector<ScalarType>;
        struct DiffRecord {
            size_t idFirstOperand;
            ExpressionType source;
        };

        Utils::Array<DiffRecord> records;
        Utils::Array<DiffScalar> operands;
        VectorType values;
        VectorType tangents;
    public:
        TraceSegment();
        TraceSegment(size_t size);
        TraceSegment(const TraceSegment&) = default;
        TraceSegment(TraceSegment&&) noexcept = default;
        ~TraceSegment() = default;
        /* Operators */
        TraceSegment& operator=(TraceSegment obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void reverse();
        void reverse(DiffScalar from);
        void reverse(DiffScalar from, DiffScalar to);
        void zero_grad(DiffScalar from);
        void zero_grad(DiffScalar from, DiffScalar to);
        void forget(DiffScalar from);
        void squeeze();
        void swap(TraceSegment& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return records.getLength(); }
        [[nodiscard]] size_t getCapacity() const noexcept { return records.getCapacity(); }
        [[nodiscard]] bool empty() const noexcept { return records.empty(); }
        [[nodiscard]] bool full() const noexcept { return records.full(); }
        [[nodiscard]] bool isFound(DiffScalar s) const noexcept { return find(s) < getLength(); }
        /* Static members */
        [[nodiscard]] constexpr static unsigned int numOperand(ExpressionType type);
    private:
        /* Operations */
        void reverse(size_t fromIndex, size_t toIndex);
        template<size_t Size>
        void reverseMulAdd(DiffScalar firstOpX, DiffScalar firstOpY, DiffScalar firstOpZ, size_t traceId);
        /* Getters */
        [[nodiscard]] size_t find(DiffScalar s) const noexcept;
    };

    template<class ScalarType>
    TraceSegment<ScalarType>::TraceSegment() : TraceSegment(DefaultSize) {}

    template<class ScalarType>
    TraceSegment<ScalarType>::TraceSegment(size_t size) {
        assert(size >= DefaultSize && "[Error]: Allocate a small segment maybe bad to performance");
        records.reserve(size);
        operands.reserve(3 * size); //MulAdd operation is 3-operand
        values.reserve(size);
        tangents.reserve(size);
    }

    template<class ScalarType>
    inline void TraceSegment<ScalarType>::reverse() {
        assert(!empty());
        reverse(getLength() - 1, 0);
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::reverse(DiffScalar from) {
        assert(!empty());
        size_t fromIndex = find(from);
        const bool isFromNotFound = fromIndex >= getLength();
        if (isFromNotFound)
            fromIndex = getLength() - 1;
        reverse(fromIndex, 0);
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::reverse(DiffScalar from, DiffScalar to) {
        assert(!empty());
        size_t fromIndex = find(from);
        const bool isFromNotFound = fromIndex >= getLength();
        if (isFromNotFound)
            fromIndex = getLength() - 1;

        size_t toIndex = find(to);
        const bool isToNotFound = toIndex >= getLength();
        if (isToNotFound)
            toIndex = 0;

        reverse(fromIndex, toIndex);
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::zero_grad(DiffScalar from) {
        size_t fromIndex = find(from);
        const bool isFromNotFound = fromIndex >= getLength();
        if (isFromNotFound)
            fromIndex = 0;
        auto segment = tangents.segment(fromIndex, getLength());
        segment = ScalarType(0);
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::zero_grad(DiffScalar from, DiffScalar to) {
        size_t fromIndex = find(from);
        const bool isFromNotFound = fromIndex >= getLength();
        if (isFromNotFound)
            fromIndex = 0;

        size_t toIndex = find(to);
        const bool isToNotFound = toIndex >= getLength();
        if (isToNotFound)
            toIndex = getLength();
        assert(fromIndex < toIndex && "[Error]: Invalid range");
        auto segment = tangents.segment(fromIndex, toIndex);
        segment = ScalarType(0);
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::forget(DiffScalar from) {
        assert(isFound(from) && "[Error]: forgeting a non existent, this may be a bug");
        const size_t index = find(from);
        operands.resize(records[index].idFirstOperand);
        records.resize(index);
        values.resize(index);
        tangents.resize(index);
    }
    /**
     * FIXME: squeeze values and tangents may invalidate pointers, so the unused memory cannot be freed
     */
    template<class ScalarType>
    void TraceSegment<ScalarType>::squeeze() {
        records.squeeze();
        operands.squeeze();
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::swap(TraceSegment& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        records.swap(obj.records);
        operands.swap(obj.operands);
        values.swap(obj.values);
        tangents.swap(obj.tangents);
    }

    template<class ScalarType>
    constexpr unsigned int TraceSegment<ScalarType>::numOperand(ExpressionType type) {
        switch (type) {
            case ExpressionType::Set: return 0;
            case ExpressionType::Assign: return 1;
            case ExpressionType::Minus: return 1;
            case ExpressionType::Add: return 2;
            case ExpressionType::Sub: return 2;
            case ExpressionType::Mul: return 2;
            case ExpressionType::Div: return 2;
            case ExpressionType::MulAdd2: [[fallthrough]];
            case ExpressionType::MulAdd4: [[fallthrough]];
            case ExpressionType::MulAdd8: return 3;
            case ExpressionType::More: return 2;
            case ExpressionType::MoreEq: return 2;
            case ExpressionType::Reciprocal: return 1;
            case ExpressionType::Sqrt: return 1;
            case ExpressionType::Cbrt: return 1;
            case ExpressionType::Abs: return 1;
            case ExpressionType::Relu: return 1;
            case ExpressionType::Square: return 1;
            case ExpressionType::Ln: return 1;
            case ExpressionType::Exp: return 1;
            case ExpressionType::Pow: return 2;
            case ExpressionType::Sin: return 1;
            case ExpressionType::Cos: return 1;
            case ExpressionType::ArcCos: return 1;
            default: [[unlikely]]
                throw std::invalid_argument("[Error]: Unrecognized type");
        }
    }

    template<class ScalarType>
    void TraceSegment<ScalarType>::reverse(size_t fromIndex, size_t toIndex) {
        assert(toIndex <= fromIndex && "[Error]: Unexpected index");
        for (size_t i = fromIndex; i >= toIndex && i <= fromIndex; --i) {
            const DiffRecord record = records[i];
            const ScalarType& tangent = tangents[i];
            if (tangent.isZero() || record.source == ExpressionType::Set)
                continue;
            
            const ScalarType& value = values[i];
            const size_t idFirstOperand = record.idFirstOperand;
            const ExpressionType source = record.source;
            DiffScalar operandX = operands[idFirstOperand];
            ScalarType& tangentX = operandX.getTangent();
            /* Unitary Operations */ {
                switch (source) {
                    case ExpressionType::Assign:
                    case ExpressionType::Minus:
                        tangentX += tangent * ScalarType(source == ExpressionType::Assign ? 1.0 : -1.0);
                        continue;
                    case ExpressionType::Reciprocal:
                        tangentX -= tangent * square(value);
                        continue;
                    case ExpressionType::Sqrt:
                        tangentX += tangent / value * ScalarType(0.5);
                        continue;
                    case ExpressionType::Cbrt:
                        tangentX += tangent / (square(value) * ScalarType(3));
                        continue;
                    case ExpressionType::Abs:
                        tangentX += operandX.getValue().isPositive() ? tangent : -tangent;
                        continue;
                    case ExpressionType::Relu:
                        tangentX += operandX.getValue().isPositive() ? tangent : ScalarType(0);
                        continue;
                    case ExpressionType::Square:
                        tangentX += tangent * operandX.getValue() * ScalarType(2);
                        continue;
                    case ExpressionType::Ln:
                        tangentX += tangent / operandX.getValue();
                        continue;
                    case ExpressionType::Exp:
                        tangentX += tangent * value;
                        continue;
                    case ExpressionType::Sin:
                        tangentX += tangent * cos(operandX.getValue());
                        continue;
                    case ExpressionType::Cos:
                        tangentX -= tangent * sin(operandX.getValue());
                        continue;
                    case ExpressionType::ArcCos:
                        tangentX -= tangent / sqrt(ScalarType(1) - square(operandX.getValue()));
                        continue;
                    default:;
                }
            }
            /* Binary Operations */ {
                DiffScalar operandY = operands[idFirstOperand + 1];
                ScalarType& tangentY = operandY.getTangent();
                switch (source) {
                    case ExpressionType::Add:
                    case ExpressionType::Sub:
                        tangentX += tangent;
                        tangentY += tangent * ScalarType(source == ExpressionType::Add ? 1.0 : -1.0);
                        continue;
                    case ExpressionType::Mul:
                        tangentX += tangent * operandY.getValue();
                        tangentY += tangent * operandX.getValue();
                        continue;
                    case ExpressionType::Div: {
                        const ScalarType dx = tangent * reciprocal(operandY.getValue());
                        tangentX += dx;
                        tangentY -= dx * value;
                        continue;
                    }
                    case ExpressionType::MulAdd2:
                        if constexpr (ScalarType::Option == Double) {
                            DiffScalar operandZ = operands[idFirstOperand + 2];
                            i -= 1;
                            [[maybe_unused]] const bool isInRange = i >= toIndex && i <= fromIndex;
                            assert(isInRange && "[Error]: Unexpected id");
                            reverseMulAdd<2>(operandX, operandY, operandZ, i);
                        }
                        continue;
                    case ExpressionType::MulAdd4: {
                        DiffScalar operandZ = operands[idFirstOperand + 2];
                        i -= 3;
                        [[maybe_unused]] const bool isInRange = i >= toIndex && i <= fromIndex;
                        assert(isInRange && "[Error]: Unexpected id");
                        reverseMulAdd<4>(operandX, operandY, operandZ, i);
                        continue;
                    }
                    case ExpressionType::MulAdd8: {
                        DiffScalar operandZ = operands[idFirstOperand + 2];
                        i -= 7;
                        [[maybe_unused]] const bool isInRange = i >= toIndex && i <= fromIndex;
                        assert(isInRange && "[Error]: Unexpected id");
                        reverseMulAdd<8>(operandX, operandY, operandZ, i);
                        continue;
                    }
                    default:;
                }
            }
            assert(false && "[Error]: Undefined operator for back propagation");
        }
    }

    template<class ScalarType>
    template<size_t Size>
    void TraceSegment<ScalarType>::reverseMulAdd(DiffScalar firstOpX, DiffScalar firstOpY, DiffScalar firstOpZ, size_t traceId) {
        using PacketType = SIMD<ScalarType, Size>;
        const PacketType tangent = tangents.template packet<PacketType>(traceId);
        PacketType tangentZ{};
        tangentZ.load(firstOpZ.tangent_ptr());
        tangentZ += tangent;
        tangentZ.store(firstOpZ.tangent_ptr());


        PacketType valueX{}, valueY{}, tangentX{}, tangentY{};
        valueX.load(firstOpX.value_ptr());
        valueY.load(firstOpY.value_ptr());
        tangentX.load(firstOpX.tangent_ptr());
        tangentY.load(firstOpY.tangent_ptr());
        tangentX += tangent * valueY;
        tangentY += tangent * valueX;
        valueX.store(firstOpX.value_ptr());
        valueY.store(firstOpY.value_ptr());
        tangentX.store(firstOpX.tangent_ptr());
        tangentY.store(firstOpY.tangent_ptr());
    }

    template<class ScalarType>
    size_t TraceSegment<ScalarType>::find(DiffScalar s) const noexcept {
        const auto* pValue = values.data();
        const auto* pValue1 = s.value_ptr();
        if (pValue1 < pValue)
            return values.getLength();
        const size_t index = pValue1 - pValue;
        [[maybe_unused]] const bool isValueFound = index < getLength();
        [[maybe_unused]] const bool isTangentFound = s.tangent_ptr() == tangents.data() + index;
        assert((!isValueFound || (isValueFound && isTangentFound)) && "[Error]: Bad DiffScalar");
        return index;
    }
}
