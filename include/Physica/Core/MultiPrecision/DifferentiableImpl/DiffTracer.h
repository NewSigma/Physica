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
#include "Physica/Core/Exception/NotImplementedException.h"

namespace Physica::Core {
    template<class ScalarType>
    class DiffTracer {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        struct DiffRecord {
            size_t startOperandId;
            ScalarType value;
            ScalarType tangent;
            ExpressionType source;
        };
        using RecordArray = Utils::Array<DiffRecord>;

        RecordArray records;
        Utils::Array<size_t> operands;
    public:
        ~DiffTracer() = default;
        /* Operations */
        inline DiffTracer& pushOperand(size_t operand);
        [[nodiscard]] inline size_t pushOperation(ScalarType value, ExpressionType source);
        [[nodiscard]] size_t pushOperation(ScalarType value, ScalarType tangent, ExpressionType source);
        void reverse(size_t index);
        /* Getters */
        [[nodiscard]] const RecordArray& getRecords() const noexcept { return records; }
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
    };

    template<class ScalarType>
    inline DiffTracer<ScalarType>& DiffTracer<ScalarType>::pushOperand(size_t operand) {
        assert(!records.empty() && "[Error]: Push operand to empty operation is not allowed");
        assert(operands.getLength() < (*records.crbegin()).startOperandId + numOperand((*records.crbegin()).source)
                && "[Error]: Not enough operand slot for new operand");
        operands.append(operand);
        return *this;
    }

    template<class ScalarType>
    inline size_t DiffTracer<ScalarType>::pushOperation(ScalarType value, ExpressionType source) {
        return pushOperation(value, 0, source);
    }

    template<class ScalarType>
    size_t DiffTracer<ScalarType>::pushOperation(ScalarType value, ScalarType tangent, ExpressionType source) {
        if (!records.empty()) {
            assert(operands.getLength() == (*records.crbegin()).startOperandId + numOperand((*records.crbegin()).source)
                    && "[Error]: New record cannot begin unless last record is done");
        }
        const size_t index = records.getLength();
        records.append(DiffRecord{operands.getLength(), std::move(value), tangent, source});
        return index;
    }

    template<class ScalarType>
    void DiffTracer<ScalarType>::reverse(size_t index) {
        assert(index < getNumRecord());
        assert(records[index].tangent.isZero() && "[Error]: Last reverse result may be not cleared");
        assert(records[0].source == ExpressionType::Set && "[Error]: Unexpected source");
        records[index].tangent = ScalarType(1);
        for (size_t i = index; i != 0; --i) {
            const auto& record = records[i];
            const size_t startOperandId = record.startOperandId;
            const ScalarType& value = record.value;
            const ScalarType& tangent = record.tangent;
            switch (record.source) {
                case ExpressionType::Set:
                    break;
                case ExpressionType::Assign:
                    break;
                case ExpressionType::Minus:
                    records[operands[startOperandId]].tangent += -tangent;
                    break;
                case ExpressionType::Add:
                    records[operands[startOperandId]].tangent += tangent;
                    records[operands[startOperandId + 1]].tangent += tangent;
                    break;
                case ExpressionType::Sub:
                    records[operands[startOperandId]].tangent += tangent;
                    records[operands[startOperandId + 1]].tangent -= tangent;
                    break;
                case ExpressionType::Mul: {
                    auto& x = records[operands[startOperandId]];
                    auto& y = records[operands[startOperandId + 1]];
                    x.tangent += tangent * y.value;
                    y.tangent += tangent * x.value;
                    break;
                }
                case ExpressionType::Div: {
                    auto& x = records[operands[startOperandId]];
                    auto& y = records[operands[startOperandId + 1]];
                    const ScalarType dx = tangent * reciprocal(y.value);
                    x.tangent += dx;
                    y.tangent -= dx * value;
                    break;
                }
                case ExpressionType::Reciprocal:
                    records[operands[startOperandId]].tangent -= square(value);
                    break;
                case ExpressionType::Sqrt: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent += tangent / value * ScalarType(0.5);
                    break;
                }
                case ExpressionType::Cbrt:
                    records[operands[startOperandId]].tangent += tangent / (square(value) * ScalarType(3));
                    break;
                case ExpressionType::Abs: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent += x.value.isPositive() ? tangent : -tangent;
                    break;
                }
                case ExpressionType::Square: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent += tangent * x.value * ScalarType(2);
                    break;
                }
                case ExpressionType::Ln: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent += tangent / x.value;
                    break;
                }
                case ExpressionType::Exp: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent += tangent * value;
                    break;
                }
                case ExpressionType::Sin: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent -= tangent * sin(x.value);
                    break;
                }
                case ExpressionType::Cos: {
                    auto& x = records[operands[startOperandId]];
                    x.tangent += tangent * cos(x.value);
                    break;
                }
                default: [[unlikely]]
                    throw NotImplementedException("[Error]: Undefined operator for back propagation");
            }
        }
    }

    template<class ScalarType>
    DiffTracer<ScalarType>& DiffTracer<ScalarType>::getInstance() noexcept {
        thread_local static DiffTracer<ScalarType> instance{};
        return instance;
    }
}
