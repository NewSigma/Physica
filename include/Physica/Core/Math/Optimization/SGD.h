/*
 * Copyright 2023 WeiBo He.
 *
 * This file is part of Physica.

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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"
#include "Physica/Core/MultiPrecision/Differentiable.h"

namespace Physica::Core {
    /**
     * \class SGD Stochastic gradient descent using auto diff
     */
    template<class ScalarType>
    class SGD {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");

        ScalarType learnRate;
        ScalarType meanLearnRate;
        unsigned int batchSize;
    public:
        SGD(ScalarType learnRate_, unsigned int batchSize_);
        SGD(const SGD&) = default;
        SGD(SGD&&) noexcept = default;
        ~SGD() = default;
        /* Operators */
        SGD& operator=(SGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void step(Differentiable<ScalarType, DiffMode::Reverse>& target) const;
        template<class VectorType> void step(LValueVector<VectorType>& target) const;
        template<class MatrixType> void step(LValueMatrix<MatrixType>& target) const;
        void swap(SGD& obj) noexcept;
        /* Getters */
        [[nodiscard]] const ScalarType& getLearnRate() const noexcept { return learnRate; }
        [[nodiscard]] unsigned int getBatchSize() const noexcept { return batchSize; }
        /* Setters */
        void setLearnRate(ScalarType lr);
    };

    template<class ScalarType>
    SGD<ScalarType>::SGD(ScalarType learnRate_, unsigned int batchSize_) : batchSize(batchSize_) {
        setLearnRate(learnRate_);
    }

    template<class ScalarType>
    inline void SGD<ScalarType>::step(Differentiable<ScalarType, DiffMode::Reverse>& target) const {
        target.setValue(target.getValue() - meanLearnRate * target.getTangent());
    }

    template<class ScalarType>
    template<class VectorType>
    void SGD<ScalarType>::step(LValueVector<VectorType>& target) const {
        using OtherScalar = typename VectorType::ScalarType;
        static_assert(OtherScalar::isDifferentiable, "[Error]: Vector must be differentiable");
        for (size_t i = 0; i < target.getLength(); ++i)
            step(target[i]);
    }
    
    template<class ScalarType>
    template<class MatrixType>
    void SGD<ScalarType>::step(LValueMatrix<MatrixType>& target) const {
        using OtherScalar = typename MatrixType::ScalarType;
        static_assert(OtherScalar::isDifferentiable, "[Error]: Matrix must be differentiable");
        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                step(target.refFromMajorMinor(i, j));
    }

    template<class ScalarType>
    void SGD<ScalarType>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
        meanLearnRate.swap(obj.meanLearnRate);
        std::swap(batchSize, obj.batchSize);
    }

    template<class ScalarType>
    void SGD<ScalarType>::setLearnRate(ScalarType lr) {
        learnRate = lr;
        meanLearnRate = lr / ScalarType(batchSize);
    }
}
