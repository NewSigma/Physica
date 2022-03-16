/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::Core {
    /**
     * References:
     * [1] 贾俊平, 何晓群, 金勇进. 统计学.第6版[M]. 中国人民大学出版社, 2015.239-244
     */
    template<class ScalarType>
    class VarianceAnalyzer {
        using DataSet = Utils::Array<Vector<ScalarType>>;
        DataSet data;
    public:
        VarianceAnalyzer(size_t numGroup, size_t samplePerGroup);
        ~VarianceAnalyzer() = default;
        /* Getters */
        [[nodiscard]] const DataSet& getData() const noexcept { return data; }
        [[nodiscard]] size_t getNumGroup() const noexcept { return data.getLength(); }
        [[nodiscard]] size_t getTotalNumSample() const noexcept;
        [[nodiscard]] ScalarType getParamF() const { return getMSA() / getMSE(); }
        [[nodiscard]] ScalarType relationCoeff() const;
        /* Setters */
        template<class VectorType>
        void setGroup(size_t groupId, const RValueVector<VectorType>& group) { data[groupId] = group; }
    private:
        ScalarType getSSE() const;
        ScalarType getMSE() const;
        ScalarType getSSA() const;
        ScalarType getMSA() const;
        ScalarType getTotalMean() const;
    };

    template<class ScalarType>
    VarianceAnalyzer<ScalarType>::VarianceAnalyzer(size_t numGroup, size_t samplePerGroup) : data(numGroup, samplePerGroup) {}

    template<class ScalarType>
    size_t VarianceAnalyzer<ScalarType>::getTotalNumSample() const noexcept {
        size_t count = 0;
        for (const auto& vec : data)
            count += vec.getLength();
        return count;
    }

    template<class ScalarType>
    ScalarType VarianceAnalyzer<ScalarType>::relationCoeff() const noexcept {
        const auto ssa = getSSA();
        return ssa / (getSSE() + ssa);
    }

    template<class ScalarType>
    ScalarType VarianceAnalyzer<ScalarType>::getSSE() const {
        ScalarType result = 0;
        for (const auto& vec : data) {
            const ScalarType mean_i = mean(vec);
            result += square(vec - mean_i).sum();
        }
        return result;
    }

    template<class ScalarType>
    ScalarType VarianceAnalyzer<ScalarType>::getMSE() const {
        return getSSA() / ScalarType(getTotalNumSample() - getNumGroup());
    }

    template<class ScalarType>
    ScalarType VarianceAnalyzer<ScalarType>::getMSA() const {
        const ScalarType mean_total = getTotalMean();
        ScalarType result = 0;
        for (const auto& vec : data)
            result += vec.getLength() * square(mean(vec) - mean_total);
        return result;
    }

    template<class ScalarType>
    ScalarType VarianceAnalyzer<ScalarType>::getMSA() const {
        return getSSA() / ScalarType(getNumGroup());
    }

    template<class ScalarType>
    ScalarType VarianceAnalyzer<ScalarType>::getTotalMean() const {
        ScalarType result = 0;
        for (const auto& vec : data)
            result += mean(vec) * ScalarType(vec.getLength());
        result /= getTotalNumSample();
        return result;
    }
}
