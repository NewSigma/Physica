/*
 * Copyright 2023-2026 Weibo He.
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
        VectorND<T> separates;
        T repDelta;
        T maximum;
        T minimum;
    public:
        ProbDistribution(T from, T to, size_t numBin);
        ProbDistribution(const This&) = default;
        ProbDistribution(This&&) noexcept = default;
        ~ProbDistribution() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        void operator+=(const This& pdf);
        [[nodiscard]] bool operator==(const This& other) const noexcept;
        /* Operations */
        void sample(T data);
        void sample(VectorND<T> datas);
        void clear();

        [[nodiscard]] VectorND<T> makePosition() const;
        [[nodiscard]] VectorND<T> makeDistribution() const;
        [[nodiscard]] size_t calcNumSample() const;

        const H5Group read(const H5Loc& loc, const char* name);
        H5Group write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getBucket() const noexcept { return bucket; }
        [[nodiscard]] size_t getNumBin() const noexcept { return bucket.getLength(); }
        [[nodiscard]] T getFromPoint() const noexcept { return separates.front(); }
        [[nodiscard]] T getToPoint() const noexcept { return separates.back(); }
        [[nodiscard]] T getMaximum() const noexcept { return maximum; }
        [[nodiscard]] T getMinimum() const noexcept { return minimum; }
    };

    template<Scalar T>
    ProbDistribution<T>::ProbDistribution(T from, T to, size_t numBin)
            : bucket(numBin, 0)
            , separates(VectorND<T>::linspace(from, to, numBin + 1))
            , repDelta(T(numBin) / (to - from)) {
        assert(from < to);
        clear();
    }

    template<Scalar T>
    auto ProbDistribution<T>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void ProbDistribution<T>::operator+=(const ProbDistribution<T>& pdf) {
        for (size_t i = 0; i < bucket.getLength(); ++i)
            bucket[i] += pdf.bucket[i];
    }

    template<Scalar T>
    bool ProbDistribution<T>::operator==(const ProbDistribution<T>& other) const noexcept {
        return bucket == other.bucket
            && getFromPoint() == other.getFromPoint()
            && getToPoint() == other.getToPoint()
            && maximum == other.maximum
            && minimum == other.minimum;
    }

    template<Scalar T>
    void ProbDistribution<T>::sample(T data) {
        const long index = double((data - getFromPoint()) * repDelta);
        if (data > getFromPoint() && 0 <= index && size_t(index) < getNumBin())
            bucket[index] += 1;
        maximum = std::max(data, maximum);
        minimum = std::min(data, minimum);
    }

    template<Scalar T>
    void ProbDistribution<T>::sample(VectorND<T> datas) {
        for (auto data : datas)
            sample(data);
    }

    template<Scalar T>
    void ProbDistribution<T>::clear() {
        bucket.zeros();
        maximum = std::numeric_limits<T>::lowest();
        minimum = std::numeric_limits<T>::max();
    }

    template<Scalar T>
    auto ProbDistribution<T>::makePosition() const -> VectorND<T> {
        const T delta = (getToPoint() - getFromPoint()) / T(getNumBin());
        return separates.head(getNumBin()) + (delta * 0.5);
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
    size_t ProbDistribution<T>::calcNumSample() const {
        size_t num = 0;
        for (auto elem : bucket)
            num += elem;
        return num;
    }

#ifdef PHYSICA_HDF5
    template<Scalar T>
    const H5Group ProbDistribution<T>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        T from, to;
        size_t numBin{};

        bucket.read(group, "Bucket");
        group.readAttr("NumBin", numBin);
        group.readAttr("From", from);
        group.readAttr("To", to);
        group.readAttr("Maximum", maximum);
        group.readAttr("Minimum", minimum);
        separates = VectorND<T>::linspace(from, to, numBin + 1);
        repDelta = T(numBin) / (to - from);
        return H5Group(group);
    }

    template<Scalar T>
    H5Group ProbDistribution<T>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        bucket.write(group, "Bucket");
        group.writeAttr("NumBin", bucket.getLength());
        group.writeAttr("From", getFromPoint());
        group.writeAttr("To", getToPoint());
        group.writeAttr("Maximum", maximum);
        group.writeAttr("Minimum", minimum);
        return H5Group(group);
    }
#endif

    template<Scalar T>
    void ProbDistribution<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        bucket.swap(obj.bucket);
        separates.swap(obj.separates);
        repDelta.swap(obj.repDelta);
        maximum.swap(obj.maximum);
        minimum.swap(obj.minimum);
    }
}
