/*
 * Copyright 2022-2024 Weibo He.
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
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica {
    /**
     * References:
     * [1] 贾俊平, 何晓群, 金勇进. 统计学.第6版[M]. 中国人民大学出版社, 2015.239-244
     */
    template<Scalar T>
    class VarianceAnalyzer {
        using DataSet = Array<VectorND<T>>;
        DataSet data;
    public:
        VarianceAnalyzer(size_t numGroup, size_t samplePerGroup);
        ~VarianceAnalyzer() = default;
        /* Getters */
        [[nodiscard]] const DataSet& getData() const noexcept { return data; }
        [[nodiscard]] size_t getNumGroup() const noexcept { return data.getLength(); }
        [[nodiscard]] size_t getTotalNumSample() const noexcept;
        [[nodiscard]] T getParamF() const { return getMSA() / getMSE(); }
        [[nodiscard]] T relationCoeff() const noexcept;
        /* Setters */
        void setGroup(size_t groupId, const Vector auto& group) { data[groupId] = group; }
    private:
        T getSSE() const;
        T getMSE() const;
        T getSSA() const;
        T getMSA() const;
        T getTotalMean() const;
    };

    template<Scalar T>
    VarianceAnalyzer<T>::VarianceAnalyzer(size_t numGroup, size_t samplePerGroup) : data(numGroup, samplePerGroup) {}

    template<Scalar T>
    size_t VarianceAnalyzer<T>::getTotalNumSample() const noexcept {
        size_t count = 0;
        for (const auto& vec : data)
            count += vec.getLength();
        return count;
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::relationCoeff() const noexcept {
        const auto ssa = getSSA();
        return ssa / (getSSE() + ssa);
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::getSSE() const {
        T result = 0;
        for (const auto& vec : data)
            result += square(vec - vec.mean()).sum();
        return result;
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::getMSE() const {
        return getSSA() / T(getTotalNumSample() - getNumGroup());
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::getSSA() const {
        const T mean_total = getTotalMean();
        T result = 0;
        for (const auto& vec : data)
            result += vec.getLength() * square(vec.mean() - mean_total);
        return result;
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::getMSA() const {
        return getSSA() / T(getNumGroup());
    }

    template<Scalar T>
    T VarianceAnalyzer<T>::getTotalMean() const {
        T result = 0;
        for (const auto& vec : data)
            result += vec.mean() * T(vec.getLength());
        result /= getTotalNumSample();
        return result;
    }
}
