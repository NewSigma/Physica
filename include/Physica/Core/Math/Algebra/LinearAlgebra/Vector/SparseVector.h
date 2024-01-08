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

#include "SparseVectorImpl/RSparseVector.h"

namespace Physica::Core {
    template<class ScalarType> class SparseVector;

    namespace Internal {
        template<class T>
        class Traits<SparseVector<T>> {
        public:
            using ScalarType = T;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class ScalarType>
    class SparseVector : public RSparseVector<SparseVector<ScalarType>> {
        using This = SparseVector<ScalarType>;
        using Base = RSparseVector<This>;
        using IndexArray = Utils::Array<size_t>;
        using ElemArray = Utils::Array<ScalarType>;
        using typename Base::NonZeroPair;

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
        /* Operations */
        template<class OtherDerived>
        void assignTo_sparse(LValueVector<OtherDerived>& v) const;
        void insert(size_t index, ScalarType value);
        void resize(size_t newLength);
        void reserve(size_t numNonZero);
        void swap(SparseVector& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] NonZeroPair calcNonZero(size_t index) const { return std::make_pair(indexes[index], elems[index]); }
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elems.getLength(); }
    };

    template<class ScalarType>
    SparseVector<ScalarType>::SparseVector() : length(0) {}

    template<class ScalarType>
    SparseVector<ScalarType>::SparseVector(size_t length_) : length(length_) {}

    template<class ScalarType>
    SparseVector<ScalarType>::SparseVector(size_t length_, size_t numNonZero) : SparseVector(length_) {
        assert(numNonZero < length && "[Error]: Too many non zero elements");
        reserve(numNonZero);
    }

    template<class ScalarType>
    template<class OtherDerived>
    void SparseVector<ScalarType>::assignTo_sparse(LValueVector<OtherDerived>& v) const {
        for (size_t i = 0; i < getNumNonZero(); ++i)
            v[indexes[i]] = elems[i];
    }

    template<class ScalarType>
    void SparseVector<ScalarType>::insert(size_t index, ScalarType value) {
        assert(index < getLength() && "[Error]: Index out of range");
        for (size_t i = 0; i < indexes.getLength(); ++i) {
            if (index == indexes[i]) {
                elems[i] = value;
                return;
            }
        }
        indexes.append(index);
        elems.append(std::move(value));
    }

    template<class ScalarType>
    void SparseVector<ScalarType>::resize(size_t newLength) {
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

    template<class ScalarType>
    void SparseVector<ScalarType>::reserve(size_t numNonZero) {
        indexes.reserve(numNonZero);
        elems.reserve(numNonZero);
    }

    template<class ScalarType>
    void SparseVector<ScalarType>::swap(SparseVector& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(length, obj.length);
        indexes.swap(obj.indexes);
        elems.swap(obj.elems);
    }

    template<class ScalarType>
    ScalarType SparseVector<ScalarType>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        for (size_t i = 0; i < indexes.getLength(); ++i) {
            if (index == indexes[i])
                return elems[i];
        }
        return ScalarType(0);
    }
}
