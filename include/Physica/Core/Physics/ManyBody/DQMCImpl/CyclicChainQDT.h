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
        using MatrixND = HubbardParams<T>::MatrixND;
        constexpr static int Option = MatrixOption::Col | MatrixOption::Element;

        Array2D<QDTDecomp<T>, Option> decomps;
        Array2D<bool, Option> readys;
    public:
        CyclicChainQDT() = default;
        explicit CyclicChainQDT(size_t length);
        CyclicChainQDT(const This&) = default;
        CyclicChainQDT(This&&) noexcept = default;
        ~CyclicChainQDT() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] QDTDecomp<T>& operator[](size_t i) { return decomps(i, i); }
        /* Operations */
        [[nodiscard]] const QDTDecomp<T>& multiply(size_t from, size_t to) noexcept;
        void invalidate(size_t i) noexcept;
        void invalidates() noexcept;

        void resize(size_t length);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return decomps.getRow(); }
    };

    template<Scalar T>
    CyclicChainQDT<T>::CyclicChainQDT(size_t length) {
        resize(length);
    }
    /**
     * Closed interval: [from, to]
     */
    template<Scalar T>
    auto CyclicChainQDT<T>::multiply(size_t from, size_t to) noexcept -> const QDTDecomp<T>& {
        assert(from < getLength() && to < getLength());
        auto& result = decomps(from, to);
        if (from == to || readys(from, to))
            return result;

        readys(from, to) = true;
        if (from < to) {
            for (size_t i = from; i < to; ++i) {
                if (readys(from, i) && readys(i + 1, to)) {
                    result = decomps(from, i) * decomps(i + 1, to);
                    return result;
                }
            }
            const size_t p = (from + to) / 2;
            result = multiply(from, p) * multiply(p + 1, to);
            return result;
        }

        assert(from > to);
        for (size_t i = from; i < getLength(); ++i) {
            size_t i1 = (i + 1) % getLength();
            if (readys(from, i) && readys(i1, to)) {
                result = decomps(from, i) * decomps(i1, to);
                return result;
            }
        }

        for (size_t i = 0; i < to; ++i) {
            if (readys(from, i) && readys(i + 1, to)) {
                result = decomps(from, i) * decomps(i + 1, to);
                return result;
            }
        }

        result = multiply(from, getLength() - 1) * multiply(0, to);
        return result;
    }

    template<Scalar T>
    void CyclicChainQDT<T>::invalidate(size_t i) noexcept {
        assert(i < getLength());
        for (size_t from = 0; from < getLength(); ++from) {
            for (size_t to = 0; to < getLength(); ++to) {
                if (from == to)
                    continue;

                bool cond1 = from <= i;
                bool cond2 = i <= to;
                if (from < to)
                    readys(from, to) &= !(cond1 && cond2);
                else
                    readys(from, to) &= !(cond1 || cond2);
            }
        }
    }

    template<Scalar T>
    void CyclicChainQDT<T>::invalidates() noexcept {
        readys.zeros();
    }

    template<Scalar T>
    void CyclicChainQDT<T>::resize(size_t length) {
        decomps.resize(length, length);
        readys.resize(length, length);
    }

    template<Scalar T>
    void CyclicChainQDT<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        decomps.swap(obj.decomps);
        readys.swap(obj.readys);
    }
}
