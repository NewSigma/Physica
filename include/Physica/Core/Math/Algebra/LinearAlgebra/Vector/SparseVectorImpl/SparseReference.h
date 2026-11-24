/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

namespace Physica {
    template<Scalar T>
    class SparseReference {
        using This = SparseReference<T>;
        using VectorType = SparseVector<T>;
        using IndexArray = VectorType::IndexArray;
        using ElemArray = VectorType::ElemArray;

        VectorType& v;
        size_t index;
    public:
        SparseReference(VectorType& v_, size_t index_);
        SparseReference(const This&) = delete;
        SparseReference(This&&) noexcept = delete;
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
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
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
