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

namespace Physica::Core {
    template<class ScalarType>
    class SparseReference {
        using This = SparseReference<ScalarType>;
        using VectorType = SparseVector<ScalarType>;
        using IndexArray = typename VectorType::IndexArray;
        using ElemArray = typename VectorType::ElemArray;

        VectorType& v;
        size_t index;
    public:
        SparseReference(VectorType& v_, size_t index_);
        SparseReference(const This&) = delete;
        SparseReference(This&&) = delete;
        ~SparseReference() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        void operator=(ScalarType value);
        void operator+=(ScalarType value);
    private:
        /* Operations */
        [[nodiscard]] ScalarType* findValue();
        void append(ScalarType value);
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] IndexArray& getIndexes() noexcept { return v.indexes; }
        [[nodiscard]] ElemArray& getElems() noexcept { return v.elems; }
    };

    template<class ScalarType>
    SparseReference<ScalarType>::SparseReference(VectorType& v_, size_t index_) : v(v_), index(index_) {
        assert(index < getLength() && "[Error]: Index out of range");
    }

    template<class ScalarType>
    void SparseReference<ScalarType>::operator=(ScalarType value) {
        auto* pValue = findValue();
        if (pValue != nullptr) {
            *pValue = value;
            return;
        }
        append(value);
    }

    template<class ScalarType>
    void SparseReference<ScalarType>::operator+=(ScalarType value) {
        auto* pValue = findValue();
        if (pValue != nullptr) {
            *pValue += value;
            return;
        }
        append(value);
    }

    template<class ScalarType>
    void SparseReference<ScalarType>::append(ScalarType value) {
        getIndexes().append(index);
        getElems().append(std::move(value));
    }

    template<class ScalarType>
    ScalarType* SparseReference<ScalarType>::findValue() {
        const auto& indexes = getIndexes();
        for (size_t i = 0; i < indexes.getLength(); ++i)
            if (index == indexes[i])
                return &getElems()[i];
        return nullptr;
    }
}
