/*
 * Copyright 2025 Weibo He.
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

#include "QDTDecomp.h"

namespace Physica {
    /**
     * \class CyclicChainQDT: A cyclic chain of QDT decompositions
     *
     * Maintains a Sliding Window Aggregation (SWAG) to achieve O(1) amortized update
     * Use QR method to ensure numerical stability
     */
    template<Scalar T>
    class CyclicChainQDT {
        using This = CyclicChainQDT<T>;
        using Tr = T::RealType;
    private:
        Array<QDTDecomp<T>> leaves;
        Array<QDTDecomp<T>> prefix;
        Array<QDTDecomp<T>> suffix;
        QDTDecomp<T> buffer;
        bool rebuild = true;
    public:
        CyclicChainQDT() = default;
        explicit CyclicChainQDT(size_t numSplit);
        CyclicChainQDT(const This&) = default;
        CyclicChainQDT(This&&) noexcept = default;
        ~CyclicChainQDT() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] auto& operator[](size_t i) { return leaves[i]; }
        /* Operations */
        [[nodiscard]] const QDTDecomp<T>& multiply(size_t from, size_t to) noexcept;
        void single_flip(int site, int split, Tr factor, Tr invfac) noexcept;
        void invalidate(int split) noexcept;
        void invalidates() noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumSplit() const noexcept { return leaves.getLength(); }
    private:
        void build(size_t cursor) noexcept;
        void transfer(size_t cursor) noexcept;
    };

    template<Scalar T>
    CyclicChainQDT<T>::CyclicChainQDT(size_t numSplit)
            : leaves(numSplit) {
        prefix.reserve(numSplit);
        suffix.reserve(numSplit);
    }

    template<Scalar T>
    auto CyclicChainQDT<T>::multiply(size_t from, size_t to) noexcept -> const QDTDecomp<T>& {
        assert(from < getNumSplit() && to < getNumSplit());
        assert(from == (to + 1) % getNumSplit() && "[Error]: CyclicChainQDT only supports the full ring product");
        if (rebuild)
            build(to);

        if (suffix.empty())
            return prefix.back();
        if (prefix.empty())
            return suffix.back();
        buffer = prefix.back() * suffix.back();
        return buffer;
    }

    template<Scalar T>
    void CyclicChainQDT<T>::single_flip(int site, int split, Tr factor, Tr invfac) noexcept {
        assert(split < getNumSplit());
        assert(!suffix.empty());
        leaves[split].single_flip(site, factor, invfac);
        suffix.back().single_flip(site, factor, invfac);
    }

    template<Scalar T>
    void CyclicChainQDT<T>::invalidate(int split) noexcept {
        assert(split < getNumSplit());
        const size_t cursor = (split + 1) % getNumSplit();
        if (prefix.empty())
            transfer(cursor);

        prefix.pop_back();
        if (suffix.empty())
            suffix.push_back(leaves[cursor]);
        else
            suffix.push_back(suffix.back() * leaves[cursor]);
    }

    template<Scalar T>
    void CyclicChainQDT<T>::invalidates() noexcept {
        rebuild = true;
    }

    template<Scalar T>
    void CyclicChainQDT<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        leaves.swap(obj.leaves);
        prefix.swap(obj.prefix);
        suffix.swap(obj.suffix);
        buffer.swap(obj.buffer);
        std::swap(rebuild, obj.rebuild);
    }

    template<Scalar T>
    void CyclicChainQDT<T>::build(size_t cursor) noexcept {
        const size_t numSplit = getNumSplit();
        prefix.clear();
        suffix.clear();
        for (size_t i = 0; i < numSplit; ++i) {
            const auto& leave = leaves[(cursor + 1 + i) % numSplit];
            if (i == 0)
                suffix.push_back(leave);
            else
                suffix.push_back(suffix[i - 1] * leave);
        }
        rebuild = false;
    }
    /**
     * Sends suffix to prefix
     */
    template<Scalar T>
    void CyclicChainQDT<T>::transfer(size_t cursor) noexcept {
        const size_t size = suffix.size();
        for (size_t i = 0; i < size; ++i) {
            const auto& leave = leaves[(cursor + size - i - 1) % getNumSplit()];
            if (prefix.empty())
                prefix.push_back(leave);
            else
                prefix.push_back(leave * prefix.back());
        }
        suffix.clear();
    }
}
