/*
 * Copyright 2020-2025 Weibo He.
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

#include "Physica/Core/Utils/Container/SymmArray.h"
#include "DenseMatrix.h"

namespace Physica {
    /**
     * \class MatrixChain calculates the product of matrices efficiently and stably.
     *
     * Reference:
     * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, Clifford Stein . 算法导论(第三版)[M]. 北京: 机械工业出版社, 2013:210-215
     */
    template<Scalar T>
    class MatrixChain {
        using This = MatrixChain<T>;
        using MatrixType = DenseMatrix<T>;

        Array<MatrixType> chain;
        SymmArray<size_t> price;
        SymmArray<size_t> point;
    public:
        explicit MatrixChain(size_t length);
        MatrixChain(const This&) = default;
        MatrixChain(This&&) noexcept = default;
        ~MatrixChain() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] MatrixType& operator[](size_t i) { return chain[i]; }
        /* Operations */
        void dynamicProgram();
        MatrixType multiply() const;
        MatrixType multiply(size_t from, size_t to) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return chain.getLength(); }
    private:
        size_t deltaPrice(size_t n) const noexcept;
        bool isCheaper(size_t new_p, size_t old_p, size_t n, size_t middle) const noexcept;
    };

    template<Scalar T>
    MatrixChain<T>::MatrixChain(size_t length)
            : chain(length), price(length), point(length) {
        assert(length >= 2 && "[Error]: Invalid length");
        for (size_t i = 0; i < length; ++i)
            price(i, i) = 0;
    }

    template<Scalar T>
    void MatrixChain<T>::dynamicProgram() {
        const size_t length = getLength();
        for (size_t l = 1; l < length; ++l) { //(l + 1) is the length of sub-chain.
            const size_t m_max = length - l;
            for (size_t m = 0; m < m_max; ++m) { // m is the start index of each sub-chain.
                const size_t m_end = m + l; // End index of each sub-chain.
                const size_t middle = (m + m_end) / 2;
                size_t& price_m = price(m, m_end);
                size_t& point_m = point(m, m_end);
                // Cut (m ... m_end) into (m ... n) and (n + 1 ... m_end).
                size_t n = m;
                /* Handle n = m */ {
                    price_m = deltaPrice(n) + price(n + 1, m_end);
                    point_m = n;
                    ++n;
                }
                for (; n < m_end; ++n) {
                    size_t p = price(m, n) + deltaPrice(n) + price(n + 1, m_end);
                    if (isCheaper(p, price_m, n, middle)) {
                        price_m = p;
                        point_m = n;
                    }
                }
            }
        }
    }

    template<Scalar T>
    auto MatrixChain<T>::multiply() const -> MatrixType {
        return multiply(0, getLength() - 1);
    }
    /**
     * Closed interval: [from, to]
     */
    template<Scalar T>
    auto MatrixChain<T>::multiply(size_t from, size_t to) const -> MatrixType {
        assert(from <= to);
        assert(to < getLength());
        if (from == to)
            return chain[from];
        const size_t cut_at = point(from, to);
        auto first = multiply(from, cut_at);
        auto second = multiply(cut_at + 1, to);
        return first * second;
    }

    template<Scalar T>
    void MatrixChain<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        chain.swap(obj.chain);
        std::swap(price, obj.price);
        std::swap(point, obj.point);
    }

    template<Scalar T>
    size_t MatrixChain<T>::deltaPrice(size_t n) const noexcept {
        return chain[n].getSize() * chain[n + 1].getRow();
    }

    template<Scalar T>
    bool MatrixChain<T>::isCheaper(size_t new_p, size_t old_p, size_t n, size_t middle) const noexcept {
        bool flag1 = new_p < old_p;
        bool flag2 = (new_p == old_p) && (n <= middle); // if price is equal, use bisection method to keep stability
        return flag1 == flag2;
    }
}
