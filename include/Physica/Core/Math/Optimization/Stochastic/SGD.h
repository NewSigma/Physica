/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"
#include "Physica/Core/MultiPrecision/Differentiable.h"

namespace Physica::Core {
    /**
     * \class SGD Stochastic gradient descent using auto diff
     */
    template<class PlainScalar>
    class SGD {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
    public:
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;
    private:
        constexpr static int AnyValue = 0;

        PlainScalar learnRate;
        PlainScalar meanLearnRate;
        unsigned int batchSize;
        ScalarType from;
        ScalarType to;
    public:
        SGD(PlainScalar learnRate_, unsigned int batchSize_);
        SGD(const SGD&) = default;
        SGD(SGD&&) noexcept = default;
        ~SGD() = default;
        /* Operators */
        SGD& operator=(SGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void recordBegin();
        inline void recordEnd();
        void step();
        void swap(SGD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] PlainScalar getLearnRate() const noexcept { return learnRate; }
        [[nodiscard]] PlainScalar getMeanLearnRate() const noexcept { return meanLearnRate; }
        [[nodiscard]] unsigned int getBatchSize() const noexcept { return batchSize; }
        /* Setters */
        void setLearnRate(PlainScalar lr);
    };

    template<class PlainScalar>
    SGD<PlainScalar>::SGD(PlainScalar learnRate_, unsigned int batchSize_) : batchSize(batchSize_) {
        setLearnRate(learnRate_);
    }

    template<class PlainScalar>
    inline void SGD<PlainScalar>::recordBegin() {
        from = ScalarType(AnyValue);
    }
    
    template<class PlainScalar>
    inline void SGD<PlainScalar>::recordEnd() {
        to = ScalarType(AnyValue);
    }

    template<class PlainScalar>
    void SGD<PlainScalar>::step() {
        size_t i = 0;
        ScalarType::forNode(from, to, [this, &i](ScalarType s) {
            s.setValue(s.getValue() - meanLearnRate * s.getGrad());
            i += 1;
        });
    }

    template<class PlainScalar>
    void SGD<PlainScalar>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
        meanLearnRate.swap(obj.meanLearnRate);
        std::swap(batchSize, obj.batchSize);
        from.swap(obj.from);
        to.swap(obj.to);
    }

    template<class PlainScalar>
    void SGD<PlainScalar>::setLearnRate(PlainScalar lr) {
        learnRate = lr;
        meanLearnRate = lr / PlainScalar(batchSize);
    }
}
