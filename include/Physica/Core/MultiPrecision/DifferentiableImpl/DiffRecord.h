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

namespace Physica::Core {
    template<class ScalarType>
    class DiffRecord {
        size_t startOperandId;
        ExpressionType source;
    public:
        ScalarType value;
        ScalarType tangent;
    public:
        DiffRecord() = default;
        DiffRecord(size_t startOperandId_, ExpressionType source_, ScalarType value_, ScalarType tangent_);
        DiffRecord(const DiffRecord&) = default;
        DiffRecord(DiffRecord&&) noexcept = default;
        ~DiffRecord() = default;
        /* Operators */
        DiffRecord& operator=(DiffRecord obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(DiffRecord& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getStartOperandId() const noexcept { return startOperandId; }
        [[nodiscard]] ExpressionType getSource() const noexcept { return source; }
    };

    template<class ScalarType>
    DiffRecord<ScalarType>::DiffRecord(size_t startOperandId_, ExpressionType source_, ScalarType value_, ScalarType tangent_)
        : startOperandId(startOperandId_)
        , source(source_)
        , value(value_)
        , tangent(tangent_) {}
    
    template<class ScalarType>
    void DiffRecord<ScalarType>::swap(DiffRecord& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(startOperandId, obj.startOperandId);
        std::swap(source, obj.source);
        value.swap(obj.value);
        tangent.swap(obj.tangent);
    }
}
