/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Scalar T>
    class ProbDistribution {
        using This = ProbDistribution<T>;

        Array<size_t> bucket;
        VectorND<T> seperates;
        T repDelta;
    public:
        ProbDistribution(T from, T to, size_t numBin);
        ProbDistribution(const This&) = default;
        ProbDistribution(This&&) noexcept = default;
        ~ProbDistribution() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        void operator+=(const This& pdf);
        /* Operations */
        void sample(T data);
        inline void sample(VectorND<T> datas);
        void clear();
        [[nodiscard]] VectorND<T> makePosition() const;
        [[nodiscard]] VectorND<T> makeDistribution() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getBucket() const noexcept { return bucket; }
        [[nodiscard]] size_t getNumBin() const noexcept { return bucket.getLength(); }
        [[nodiscard]] T getFromPoint() const noexcept { return seperates[0]; }
        [[nodiscard]] T getToPoint() const noexcept { return *seperates.crbegin(); }
    private:
        size_t calcNumSample() const;
    };

    template<Scalar T>
    ProbDistribution<T>::ProbDistribution(T from, T to, size_t numBin)
            : bucket(numBin, 0)
            , seperates(VectorND<T>::linspace(from, to, numBin + 1))
            , repDelta(T(numBin) / (to - from)) {
        assert(from < to);
    }

    template<Scalar T>
    void ProbDistribution<T>::sample(T data) {
        const long index = double((data - getFromPoint()) * repDelta);
        if (data > getFromPoint() && 0 <= index && size_t(index) < getNumBin())
            bucket[index] += 1;
    }

    template<Scalar T>
    inline void ProbDistribution<T>::sample(VectorND<T> datas) {
        for (auto data : datas)
            sample(data);
    }

    template<Scalar T>
    void ProbDistribution<T>::clear() {
        for (auto& elem : bucket)
            elem = 0;
    }

    template<Scalar T>
    auto ProbDistribution<T>::makePosition() const -> VectorND<T> {
        const T delta = (getToPoint() - getFromPoint()) / T(getNumBin());
        return seperates.head(getNumBin()) + (delta * 0.5);
    }

    template<Scalar T>
    auto ProbDistribution<T>::makeDistribution() const -> VectorND<T> {
        const T factor = repDelta / T(calcNumSample());
        VectorND<T> result(getNumBin());
        for (size_t i = 0; i < result.getLength(); ++i)
            result[i] = T(bucket[i]) * factor;
        return result;
    }

    template<Scalar T>
    void ProbDistribution<T>::operator+=(const ProbDistribution<T>& pdf) {
        for (size_t i = 0; i < bucket.getLength(); ++i)
            bucket[i] += pdf.bucket[i];
    }

    template<Scalar T>
    void ProbDistribution<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        bucket.swap(obj.bucket);
        seperates.swap(obj.seperates);
        repDelta.swap(obj.repDelta);
    }

    template<Scalar T>
    size_t ProbDistribution<T>::calcNumSample() const {
        size_t num = 0;
        for (auto elem : bucket)
            num += elem;
        return num;
    }
}
