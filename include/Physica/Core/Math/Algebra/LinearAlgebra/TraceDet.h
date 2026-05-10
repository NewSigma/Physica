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

#include "Vector/DenseVector.h"

namespace Physica {
    /**
     * \class TraceDet: Random trace estimation for determinants. Applies to positive definite matrix
     *
     * Reference:
     * [1] J. Stat. Comput. 77(4), 329-348 (2007);  https://doi.org/10.1080/10629360600569279
     */
    template<Scalar T>
    class TraceDet {
        using This = TraceDet<T>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;

        VectorND<T> rand;
        VectorND<T> dot;
        VectorND<T> buffer;
        Tr normInf;

        int numSample;
        int numExpand;
    public:
        TraceDet(size_t size, int numSample, int numExpend);
        TraceDet(const This&) = default;
        TraceDet(This&&) noexcept = default;
        ~TraceDet() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R>
        [[nodiscard]] T compute(const Matrix auto& m);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return rand.getLength(); }
    private:
        void resize(size_t size);
    };

    template<Scalar T>
    TraceDet<T>::TraceDet(size_t size, int numSample, int numExpend)
            : rand(size), dot(size), buffer(size), numSample(numSample), numExpand(numExpend) {
        assert(numSample > 0);
        assert(numExpand > 0);
    }

    template<Scalar T>
    template<RNG R>
    T TraceDet<T>::compute(const Matrix auto& m) {
        assert(m.isSquare());
        if (m.getRow() != getOrder())
            resize(m.getRow());
        normInf = m.normInf();

        VectorND<T> terms(numExpand);
        for (int i = 0; i < numSample; ++i) {
            rand.template random_normal<R>();
            dot = rand;
            const Tr norm2 = rand.squaredNorm();

            for (int k = 0; k < numExpand; ++k) {
                buffer = (reciprocal(normInf) * m) * dot;
                dot -= buffer;
                terms[k].toNextMean(i, rand.conjugate() * dot / (norm2 * Trv(k + 1)));
            }
        }
        return T(getOrder()) * (ln(normInf) - terms.sum());
    }

    template<Scalar T>
    void TraceDet<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        rand.swap(obj.rand);
        normInf.swap(obj.normInf);

        std::swap(numSample, obj.numSample);
        std::swap(numExpand, obj.numExpand);
    }

    template<Scalar T>
    void TraceDet<T>::resize(size_t size) {
        rand.resize(size);
        dot.resize(size);
        buffer.resize(size);
    }
}
