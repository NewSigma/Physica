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
     * \class CyclicChainQDT: A chain of QDT decompositions
     *
     * Use dynamic programming to improve efficiency
     * Use QR method to ensure numerical stability
     */
    template<Scalar T>
    class CyclicChainQDT {
        using This = CyclicChainQDT<T>;

        using Tr = T::RealType;
        using OptionalQDT = std::optional<QDTDecomp<T>>;
    private:
        Array2D<OptionalQDT> decomps;
        Array<Array<OptionalQDT*>> chainBuffer;
    public:
        CyclicChainQDT() = default;
        explicit CyclicChainQDT(size_t numSplit);
        CyclicChainQDT(const This&) = default;
        CyclicChainQDT(This&&) noexcept = default;
        ~CyclicChainQDT() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] auto& operator[](size_t i) { return decomps(i, i); }
        /* Operations */
        [[nodiscard]] const QDTDecomp<T>& multiply(size_t from, size_t to) noexcept;
        void single_flip(int site, int split, Tr factor, Tr invfac) noexcept;
        void invalidate(int split) noexcept;
        void invalidates() noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getDecomps() const noexcept { return decomps; }
        [[nodiscard]] size_t getNumSplit() const noexcept { return decomps.getRow(); }
    private:
        [[nodiscard]] bool inRange(size_t from, size_t to, size_t split) const noexcept;
    };

    template<Scalar T>
    CyclicChainQDT<T>::CyclicChainQDT(size_t numSplit)
            : decomps(numSplit, numSplit)
            , chainBuffer(numSplit) {
        for (size_t split = 0; split < numSplit; ++split) {
            // Buffering reduces time complexity from O(n^2) to O(n).
            auto& subchains = chainBuffer[split];
            subchains.reserve(numSplit);
            for (size_t from = 0; from < numSplit; ++from) {
                for (size_t to = 0; to < numSplit; ++to) {
                    if (inRange(from, to, split))
                        subchains.append(decomps.data_ptr(from, to));
                }
            }
            subchains.squeeze();
            std::ranges::sort(subchains);
        }
    }
    /**
     * Closed interval: [from, to]
     */
    template<Scalar T>
    auto CyclicChainQDT<T>::multiply(size_t from, size_t to) noexcept -> const QDTDecomp<T>& {
        assert(from < getNumSplit() && to < getNumSplit());
        auto& result = decomps(from, to);
        if (result)
            return *result;

        if (from < to) {
            for (size_t i = from; i < to; ++i) {
                if (decomps(from, i) && decomps(i + 1, to)) {
                    result = (*decomps(from, i)) * (*decomps(i + 1, to));
                    return *result;
                }
            }
            const size_t p = (from + to) / 2;
            result = multiply(from, p) * multiply(p + 1, to);
            return *result;
        }

        assert(from > to);
        for (size_t i = from; i < getNumSplit(); ++i) {
            size_t i1 = (i + 1) % getNumSplit();
            if (decomps(from, i) && decomps(i1, to)) {
                result = (*decomps(from, i)) * (*decomps(i1, to));
                return *result;
            }
        }

        for (size_t i = 0; i < to; ++i) {
            if (decomps(from, i) && decomps(i + 1, to)) {
                result = (*decomps(from, i)) * (*decomps(i + 1, to));
                return *result;
            }
        }

        result = multiply(from, getNumSplit() - 1) * multiply(0, to);
        return *result;
    }

    template<Scalar T>
    void CyclicChainQDT<T>::single_flip(int site, int split, Tr factor, Tr invfac) noexcept {
        assert(split < getNumSplit());
        assert(scalarNear(factor * invfac, Tr(1), std::numeric_limits<T>::epsilon() * 10) && "[Error]: Invalid argument");
        for (size_t from = 0; from < getNumSplit(); ++from) {
            bool ready = decomps(from, split) || (from == split);
            if (ready)
                decomps(from, split)->single_flip(site, factor, invfac);
        }
    }

    template<Scalar T>
    void CyclicChainQDT<T>::invalidate(int split) noexcept {
        assert(split < getNumSplit());
        for (auto* optional : chainBuffer[split])
            optional->reset();
    }

    template<Scalar T>
    void CyclicChainQDT<T>::invalidates() noexcept {
        for (auto& elem : decomps.asArray())
            elem.reset();
    }

    template<Scalar T>
    void CyclicChainQDT<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        decomps.swap(obj.decomps);
        chainBuffer.swap(obj.chainBuffer);
    }
    /**
     * \returns true if split is within [from, to)
     */
    template<Scalar T>
    bool CyclicChainQDT<T>::inRange(size_t from, size_t to, size_t split) const noexcept {
        if (from == to)
            return false;
        bool cond1 = from <= split;
        bool cond2 = split < to;
        bool inRange1 = (from < to) && cond1 && cond2;
        bool inRange2 = !(from < to) && (cond1 || cond2);
        return inRange1 || inRange2;
    }
}
