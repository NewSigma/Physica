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

#include "SparseVectorImpl/RSparseVector.h"

namespace Physica::Core {
    template<Scalar T> class SparseReference;

    template<Scalar T>
    class SparseVector : public RSparseVector<SparseVector<T>> {
        using This = SparseVector<T>;
        using Base = RSparseVector<This>;
        using SparseRef = SparseReference<T>;
        using typename Base::NonZeroPair;
    public:
        using IndexArray = Array<size_t>;
        using ElemArray = Array<T>;
    private:
        size_t length;
        IndexArray indexes;
        ElemArray elems;
    public:
        SparseVector();
        SparseVector(size_t length_);
        SparseVector(size_t length_, size_t numNonZero);
        SparseVector(const SparseVector&) = default;
        SparseVector(SparseVector&&) noexcept = default;
        ~SparseVector() = default;
        /* Operators */
        SparseVector& operator=(SparseVector obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] SparseRef operator[](size_t index) { return SparseRef(*this, index); }
        /* Operations */
        template<LVector V>
        void assignTo_sparse(V& v) const;
        void resize(size_t newLength);
        void reserve(size_t numNonZero);
        inline void clear();
        void swap(SparseVector& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] NonZeroPair calcNonZero(size_t index) const { return std::make_pair(indexes[index], elems[index]); }
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elems.getLength(); }
        /* Friends */
        friend class SparseReference<T>;
    };

    template<Scalar T>
    SparseVector<T>::SparseVector() : length(0) {}

    template<Scalar T>
    SparseVector<T>::SparseVector(size_t length_) : length(length_) {}

    template<Scalar T>
    SparseVector<T>::SparseVector(size_t length_, size_t numNonZero) : SparseVector(length_) {
        assert(numNonZero <= length && "[Error]: Too many non zero elements");
        reserve(numNonZero);
    }

    template<Scalar T>
    template<LVector V>
    void SparseVector<T>::assignTo_sparse(V& v) const {
        for (size_t i = 0; i < getNumNonZero(); ++i)
            v[indexes[i]] = elems[i];
    }

    template<Scalar T>
    void SparseVector<T>::resize(size_t newLength) {
        if (newLength < length) {
            IndexArray newIndexes;
            ElemArray newElems;
            newIndexes.reserve(indexes.getLength());
            newElems.reserve(elems.getLength());
            for (size_t i = 0; i < indexes.getLength(); ++i) {
                const size_t index = indexes[i];
                if (index < newLength) {
                    newIndexes.grow(index);
                    newElems.grow(elems[i]);
                }
            }
            indexes.swap(newIndexes);
            elems.swap(newElems);
        }
        length = newLength;
    }

    template<Scalar T>
    void SparseVector<T>::reserve(size_t numNonZero) {
        indexes.reserve(numNonZero);
        elems.reserve(numNonZero);
    }

    template<Scalar T>
    inline void SparseVector<T>::clear() {
        indexes.clear();
        elems.clear();
    }

    template<Scalar T>
    void SparseVector<T>::swap(SparseVector& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(length, obj.length);
        indexes.swap(obj.indexes);
        elems.swap(obj.elems);
    }

    template<Scalar T>
    T SparseVector<T>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        for (size_t i = 0; i < indexes.getLength(); ++i) {
            if (index == indexes[i])
                return elems[i];
        }
        return T(0);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<SparseVector<T>> {
    public:
        using ScalarType = T;
        constexpr static size_t SizeAtCompile = Dynamic;
    };
}

#include "SparseVectorImpl/SparseReference.h"
