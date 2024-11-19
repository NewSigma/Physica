/*
 * Copyright 2024 Weibo He.
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
    template<Scalar T>
    class SparseReference {
        using This = SparseReference<T>;
        using VectorType = SparseVector<T>;
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
        void operator=(T value);
        void operator+=(T value);
    private:
        /* Operations */
        [[nodiscard]] T* findValue();
        void append(T value);
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] IndexArray& getIndexes() noexcept { return v.indexes; }
        [[nodiscard]] ElemArray& getElems() noexcept { return v.elems; }
    };

    template<Scalar T>
    SparseReference<T>::SparseReference(VectorType& v_, size_t index_) : v(v_), index(index_) {
        assert(index < getLength() && "[Error]: Index out of range");
    }

    template<Scalar T>
    void SparseReference<T>::operator=(T value) {
        auto* pValue = findValue();
        if (pValue != nullptr) {
            *pValue = value;
            return;
        }
        append(value);
    }

    template<Scalar T>
    void SparseReference<T>::operator+=(T value) {
        auto* pValue = findValue();
        if (pValue != nullptr) {
            *pValue += value;
            return;
        }
        append(value);
    }

    template<Scalar T>
    void SparseReference<T>::append(T value) {
        getIndexes().append(index);
        getElems().append(std::move(value));
    }

    template<Scalar T>
    T* SparseReference<T>::findValue() {
        const auto& indexes = getIndexes();
        for (size_t i = 0; i < indexes.getLength(); ++i)
            if (index == indexes[i])
                return &getElems()[i];
        return nullptr;
    }
}
