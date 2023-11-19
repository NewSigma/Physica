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

#include "Physica/Core/Exception/NotImplementedException.h"
#include "Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class ScalarType>
    class DiffTracer {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
    public:
        using DiffScalar = Differentiable<ScalarType, DiffMode::Reverse>;
        using VectorType = Vector<ScalarType>;
    private:
        struct DiffRecord {
            size_t startOperandId;
            ExpressionType source;
        };

        Utils::Array<DiffRecord> records;
        Utils::Array<DiffScalar> operands;
        VectorType values;
        VectorType tangents;
    public:
        ~DiffTracer() = default;
        /* Operations */
        inline void pushOperand(const DiffScalar& operand);
        template<size_t Size>
        inline void pushOperand(const DiffScalar (&operand)[Size]);

        [[nodiscard]] inline DiffScalar pushOperation(ScalarType value, ExpressionType source);
        [[nodiscard]] DiffScalar pushOperation(ScalarType value, ScalarType tangent, ExpressionType source);
        template<size_t Size>
        [[nodiscard]] DiffScalar pushOperation(const SIMD<ScalarType, Size>& simd, ExpressionType source);

        void reverse(size_t from, size_t to = 0);
        inline void zero_grad(size_t from);
        void zero_grad(size_t from, size_t to);
        void forget(size_t fromIndex = 0);
        void squeeze();
        void release();

        [[nodiscard]] inline bool checkLastOpDone() const;
        /* Getters */
        [[nodiscard]] VectorType& getValues() noexcept { return values; }
        [[nodiscard]] const VectorType& getValues() const noexcept { return values; }
        [[nodiscard]] VectorType& getTangents() noexcept { return tangents; }
        [[nodiscard]] const VectorType& getTangents() const noexcept { return tangents; }
        [[nodiscard]] const Utils::Array<DiffRecord>& getRecords() const noexcept { return records; }
        [[nodiscard]] size_t getNumRecord() const noexcept { return records.getLength(); }
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
        template<size_t Size>
        void checkOperands(const DiffScalar (&operand)[Size]);
        /* Static members */
        [[nodiscard]] constexpr static unsigned int numOperand(ExpressionType type);
    };

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::pushOperand(const DiffScalar& operand) {
        assert(!records.empty() && "[Error]: Push operand to empty operation is not allowed");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");
        assert(operand.getTraceIndex() < records.getLength() && "[Error]: This operand is not registered");
        operands.append(operand);
    }

    template<class ScalarType>
    template<size_t Size>
    inline void DiffTracer<ScalarType>::pushOperand(const DiffScalar (&operand)[Size]) {
        assert(!records.empty() && "[Error]: Push operand to empty operation is not allowed");
        assert(!checkLastOpDone() && "[Error]: Not enough operand slot for this operand");
        checkOperands(operand);

        const size_t length = operands.getLength();
        const size_t newLength = length + Size;
        operands.reserve(newLength);
        operands.setLength(newLength);

        static_assert(sizeof(DiffScalar) == 16, "[Error]: Unexpected size, the optimization is not available");
        using DoubleType = Scalar<Double>;
        auto* const p = operands.data() + length;
        if constexpr (Size % 4 == 0 && Internal::Instrset::hasAVX512()) {
            for (size_t i = 0; i < Size; i += 4) {
                SIMD<DoubleType, 8> temp{};
                temp.load(reinterpret_cast<const DoubleType*>(operand + i));
                temp.store(reinterpret_cast<DoubleType*>(p + i));
            }
        }
        else if constexpr (Size % 2 == 0 && Internal::Instrset::hasAVX()) {
            for (size_t i = 0; i < Size; i += 2) {
                SIMD<DoubleType, 4> temp{};
                temp.load(reinterpret_cast<const DoubleType*>(operand + i));
                temp.store(reinterpret_cast<DoubleType*>(p + i));
            }
        }
        else {
            for (size_t i = 0; i < Size; ++i)
                p[i] = operand[i];
        }
    }

    template<class ScalarType>
    inline typename DiffTracer<ScalarType>::DiffScalar
    DiffTracer<ScalarType>::pushOperation(ScalarType value, ExpressionType source) {
        return pushOperation(value, 0, source);
    }

    template<class ScalarType>
    typename DiffTracer<ScalarType>::DiffScalar
    DiffTracer<ScalarType>::pushOperation(ScalarType value, ScalarType tangent, ExpressionType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        const size_t index = getNumRecord();
        if (records.full()) {
            records.doubleSpace();
            values.doubleSpace();
            tangents.doubleSpace();
        }
        records.grow(DiffRecord{operands.getLength(), source});
        values.grow(std::move(value));
        tangents.grow(std::move(tangent));
        return DiffScalar(*this, index);
    }

    template<class ScalarType>
    template<size_t Size>
    typename DiffTracer<ScalarType>::DiffScalar
    DiffTracer<ScalarType>::pushOperation(const SIMD<ScalarType, Size>& simd, ExpressionType source) {
        assert(checkLastOpDone() && "[Error]: New record cannot begin unless last record is done");
        const size_t oldNumRecord = getNumRecord();
        const size_t newNumRecord = oldNumRecord + Size;
        {
            records.reserve(newNumRecord);
            records.setLength(newNumRecord);
            const size_t length = operands.getLength();
            const unsigned int numOp = numOperand(source);
            for (size_t i = 0; i < Size; ++i)
                records[oldNumRecord + i] = DiffRecord{length + i * numOp, source};
        }
        values.reserve(newNumRecord);
        simd.store(values.data_ptr(oldNumRecord));
        values.setLength(newNumRecord);

        SIMD<ScalarType, Size> empty(0);
        tangents.reserve(newNumRecord);
        empty.store(tangents.data_ptr(oldNumRecord));
        tangents.setLength(newNumRecord);
        return DiffScalar(*this, oldNumRecord);
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse(size_t from, size_t to) {
        assert(from < getNumRecord() && "[Error]: Index overflow");
        assert(to <= from);
        for (size_t i = from; i >= to && i <= from; --i) {
            const DiffRecord record = records[i];
            const ScalarType& tangent = tangents[i];
            if (tangent.isZero() || record.source == ExpressionType::Set)
                continue;
            
            const ScalarType& value = values[i];
            const size_t startOperandId = record.startOperandId;
            const ExpressionType source = record.source;
            DiffScalar operandX = operands[startOperandId];
            auto& tracerX = operandX.getTracer();
            const size_t indexX = operandX.getTraceIndex();
            ScalarType& tangentX = tracerX.tangents[indexX];
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
                        tangentX += tracerX.values[indexX].isPositive() ? tangent : -tangent;
                        continue;
                    case ExpressionType::Relu:
                        tangentX += tracerX.values[indexX].isPositive() ? tangent : ScalarType(0);
                        continue;
                    case ExpressionType::Square:
                        tangentX += tangent * tracerX.values[indexX] * ScalarType(2);
                        continue;
                    case ExpressionType::Ln:
                        tangentX += tangent / tracerX.values[indexX];
                        continue;
                    case ExpressionType::Exp:
                        tangentX += tangent * value;
                        continue;
                    case ExpressionType::Sin:
                        tangentX += tangent * cos(tracerX.values[indexX]);
                        continue;
                    case ExpressionType::Cos:
                        tangentX -= tangent * sin(tracerX.values[indexX]);
                        continue;
                    default:;
                }
            }
            /* Binary Operations */ {
                DiffScalar operandY = operands[startOperandId + 1];
                auto& tracerY = operandY.getTracer();
                const size_t indexY = operandY.getTraceIndex();
                ScalarType& tangentY = tracerY.tangents[indexY];
                switch (source) {
                    case ExpressionType::Add:
                    case ExpressionType::Sub:
                        tangentX += tangent;
                        tangentY += tangent * ScalarType(source == ExpressionType::Add ? 1.0 : -1.0);
                        continue;
                    case ExpressionType::MulAdd: {
                        DiffScalar operandZ = operands[startOperandId + 2];
                        ScalarType& tangentZ = operandZ.getTracer().tangents[operandZ.getTraceIndex()];
                        tangentZ += tangent;
                        [[fallthrough]];
                    }
                    case ExpressionType::Mul:
                        tangentX += tangent * tracerY.values[indexY];
                        tangentY += tangent * tracerX.values[indexX];
                        continue;
                    case ExpressionType::Div: {
                        const ScalarType dx = tangent * reciprocal(tracerY.values[indexY]);
                        tangentX += dx;
                        tangentY -= dx * value;
                        continue;
                    }
                    default:;
                }
            }
            assert(false && "[Error]: Undefined operator for back propagation");
        }
    }

    template<class ScalarType>
    inline void DiffTracer<ScalarType>::zero_grad(size_t from) {
        zero_grad(from, getNumRecord());
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::zero_grad(size_t from, size_t to) {
        assert(from < to && "[Error]: Invalid argument");
        assert(to <= getNumRecord() && "[Error]: Index overflow");
        auto seg = tangents.segment(from, to);
        seg = ScalarType(0);
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::forget(size_t fromIndex) {
        assert(fromIndex < getNumRecord() && "[Error]: forgeting a non existent, this may be a bug");
        operands.resize(records[fromIndex].startOperandId);
        records.resize(fromIndex);
        values.resize(fromIndex);
        tangents.resize(fromIndex);
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::squeeze() {
        records.squeeze();
        operands.squeeze();
        values.squeeze();
        tangents.squeeze();
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::release() {
        forget();
        squeeze();
    }

    template<class ScalarType>
    inline bool DiffTracer<ScalarType>::checkLastOpDone() const {
        if (records.empty())
            return true;
        return operands.getLength() >= (*records.crbegin()).startOperandId + numOperand((*records.crbegin()).source);
    }

    template<class ScalarType>
    DiffTracer<ScalarType>& DiffTracer<ScalarType>::getInstance() noexcept {
        thread_local static DiffTracer<ScalarType> instance{};
        return instance;
    }

    template<class ScalarType>
    template<size_t Size>
    void DiffTracer<ScalarType>::checkOperands([[maybe_unused]] const DiffScalar (&operand)[Size]) {
        for (size_t i = 0; i < Size; ++i) {
            assert(operand[i].getTraceIndex() < records.getLength() && "[Error]: This operand is not registered");
        }
    }

    template<class ScalarType>
    constexpr unsigned int DiffTracer<ScalarType>::numOperand(ExpressionType type) {
        switch (type) {
            case ExpressionType::Set: return 0;
            case ExpressionType::Assign: return 1;
            case ExpressionType::Minus: return 1;
            case ExpressionType::Add: return 2;
            case ExpressionType::Sub: return 2;
            case ExpressionType::Mul: return 2;
            case ExpressionType::Div: return 2;
            case ExpressionType::MulAdd: return 3;
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
            default: [[unlikely]]
                throw std::invalid_argument("[Error]: Unrecognized type");
        }
    }
}
