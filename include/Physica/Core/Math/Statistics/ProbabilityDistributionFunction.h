/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class ScalarType>
    class ProbabilityDistributionFunction {
        using This = ProbabilityDistributionFunction;
        using BucketType = Array<size_t>;
        using VectorType = Vector<ScalarType>;

        BucketType bucket;
        VectorType seperates;
        ScalarType repDelta;
    public:
        ProbabilityDistributionFunction(ScalarType from, ScalarType to, size_t numBin);
        ProbabilityDistributionFunction(const This&) = default;
        ProbabilityDistributionFunction(This&&) noexcept = default;
        ~ProbabilityDistributionFunction() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        void operator+=(const This& pdf);
        /* Operations */
        void sample(ScalarType data);
        inline void sample(VectorType datas);
        void clear();
        [[nodiscard]] VectorType makePosition() const;
        [[nodiscard]] VectorType makeDistribution() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const BucketType& getBucket() const noexcept { return bucket; }
        [[nodiscard]] size_t getNumBin() const noexcept { return bucket.getLength(); }
        [[nodiscard]] ScalarType getFromPoint() const noexcept { return seperates[0]; }
        [[nodiscard]] ScalarType getToPoint() const noexcept { return *seperates.crbegin(); }
    private:
        size_t calcNumSample() const;
    };

    template<class ScalarType>
    ProbabilityDistributionFunction<ScalarType>::ProbabilityDistributionFunction(ScalarType from, ScalarType to, size_t numBin)
            : bucket(numBin, 0)
            , seperates(VectorType::linspace(from, to, numBin + 1))
            , repDelta(ScalarType(numBin) / (to - from)) {
        assert(from < to);
    }

    template<class ScalarType>
    void ProbabilityDistributionFunction<ScalarType>::sample(ScalarType data) {
        const long index = double((data - getFromPoint()) * repDelta);
        if (data > getFromPoint() && 0 <= index && size_t(index) < getNumBin())
            bucket[index] += 1;
    }

    template<class ScalarType>
    inline void ProbabilityDistributionFunction<ScalarType>::sample(VectorType datas) {
        for (auto data : datas)
            sample(data);
    }

    template<class ScalarType>
    void ProbabilityDistributionFunction<ScalarType>::clear() {
        for (auto& elem : bucket)
            elem = 0;
    }

    template<class ScalarType>
    typename ProbabilityDistributionFunction<ScalarType>::VectorType
    ProbabilityDistributionFunction<ScalarType>::makePosition() const {
        const ScalarType delta = (getToPoint() - getFromPoint()) / ScalarType(getNumBin());
        return seperates.head(getNumBin()) + (delta * 0.5);
    }

    template<class ScalarType>
    typename ProbabilityDistributionFunction<ScalarType>::VectorType
    ProbabilityDistributionFunction<ScalarType>::makeDistribution() const {
        const ScalarType factor = repDelta / ScalarType(calcNumSample());
        VectorType result(getNumBin());
        for (size_t i = 0; i < result.getLength(); ++i)
            result[i] = ScalarType(bucket[i]) * factor;
        return result;
    }

    template<class ScalarType>
    ProbabilityDistributionFunction<ScalarType>& ProbabilityDistributionFunction<ScalarType>::operator=(
            ProbabilityDistributionFunction<ScalarType> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    void ProbabilityDistributionFunction<ScalarType>::operator+=(const ProbabilityDistributionFunction<ScalarType>& pdf) {
        for (size_t i = 0; i < bucket.getLength(); ++i)
            bucket[i] += pdf.bucket[i];
    }

    template<class ScalarType>
    void ProbabilityDistributionFunction<ScalarType>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        bucket.swap(obj.bucket);
        seperates.swap(obj.seperates);
        repDelta.swap(obj.repDelta);
    }

    template<class ScalarType>
    size_t ProbabilityDistributionFunction<ScalarType>::calcNumSample() const {
        size_t num = 0;
        for (auto elem : bucket)
            num += elem;
        return num;
    }
}
