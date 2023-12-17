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

#include <memory>
#include "SGD.h"

namespace Physica::Core {
    template<class PlainScalar>
    class MomentumSGD : public SGD<PlainScalar> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        using Base = SGD<PlainScalar>;
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;
        constexpr static int AnyValue = 0;

        Vector<PlainScalar> lastGrad;
        PlainScalar momentum;
        ScalarType from;
        ScalarType to;
    public:
        MomentumSGD(PlainScalar momentum_, PlainScalar learnRate, unsigned int batchSize);
        MomentumSGD(const MomentumSGD&) = default;
        MomentumSGD(MomentumSGD&&) noexcept = default;
        ~MomentumSGD() = default;
        /* Operators */
        MomentumSGD& operator=(MomentumSGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void recordBegin();
        void recordEnd();
        void step();
        void swap(MomentumSGD& __restrict obj) noexcept;
    private:
        using Base::step;
    };

    template<class PlainScalar>
    MomentumSGD<PlainScalar>::MomentumSGD(PlainScalar momentum_, PlainScalar learnRate, unsigned int batchSize)
            : Base(std::move(learnRate), batchSize)
            , momentum(std::move(momentum_)) {
        assert(momentum.isPositive() && "[Error]: Invalid momentum");
    }

    template<class PlainScalar>
    void MomentumSGD<PlainScalar>::recordBegin() {
        from = ScalarType(AnyValue);
    }
    
    template<class PlainScalar>
    void MomentumSGD<PlainScalar>::recordEnd() {
        to = ScalarType(AnyValue);
        lastGrad.resize(ScalarType::distance(from, to), PlainScalar(0));
    }

    template<class PlainScalar>
    void MomentumSGD<PlainScalar>::step() {
        size_t i = 0;
        ScalarType::forNode(from, to, [this, &i](ScalarType s) {
            lastGrad[i] = momentum * lastGrad[i] + s.getGrad();
            s.setValue(s.getValue() - Base::getMeanLearnRate() * lastGrad[i]);
            i += 1;
        });
    }

    template<class PlainScalar>
    void MomentumSGD<PlainScalar>::swap(MomentumSGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        momentum.swap(obj.momentum);
        from.swap(obj.from);
        to.swap(obj.to);
    }
}
