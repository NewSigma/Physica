/*
 * Copyright 2024 Weibo He.
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

#include <Physica/Core/MultiPrecision/Diff.h>
#include "Vector.h"

namespace Physica::Core {
    template<class PlainScalar, unsigned int Order>
    class Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>
            : public RValueVector<Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>> {
        using VectorType = Vector<PlainScalar>;
        using This = Diff<VectorType, DiffMode::Reverse, Order>;
        using Base = RValueVector<This>;
        using TracerType = DiffTracer<PlainScalar, Order>;
        using SegmentType = TraceSegment<PlainScalar, Order>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Diff();
        Diff(size_t length);
        Diff(VectorType values);
        Diff(const Diff&) = default;
        Diff(Diff&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        Diff& operator=(const Diff&) = default;
        Diff& operator=(Diff&&) noexcept = default;
        /* Operations */
        template<class RandomGenerator> inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator> inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(Diff& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return traceSeg.getLength(); }
        [[nodiscard]] const VectorType& getValue() const noexcept { return traceSeg.getValues(); }
        [[nodiscard]] const VectorType& getGrad() const noexcept { return traceSeg.getGrads(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist, RandomGenerator& gen);
        /* Friends */
        friend class device_obj<This>;
    };

    template<class PlainScalar, unsigned int Order>
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::Diff()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<class PlainScalar, unsigned int Order>
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::Diff(size_t length)
            : traceSeg(TracerType::getInstance().pushSegment(length)) {}

    template<class PlainScalar, unsigned int Order>
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::Diff(VectorType values)
            : traceSeg(TracerType::getInstance().pushSegment(std::move(values))) {}

    template<class PlainScalar, unsigned int Order>
    template<class RandomGenerator>
    inline void Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getLength(), gen);
    }

    template<class PlainScalar, unsigned int Order>
    template<class RandomGenerator>
    inline void Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getLength(), gen);
    }

    template<class PlainScalar, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline void Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getLength(), dist, gen);
    }

    template<class PlainScalar, unsigned int Order>
    inline typename Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::ScalarType
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        return traceSeg[index];
    }

    template<class PlainScalar, unsigned int Order>
    template<class RandomGenerator>
    inline Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::random_uniform(size_t len, RandomGenerator& gen) {
        return This(VectorType::random_uniform(len, gen));
    }

    template<class PlainScalar, unsigned int Order>
    template<class RandomGenerator>
    inline Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::random_normal(size_t len, RandomGenerator& gen) {
        return This(VectorType::random_normal(len, gen));
    }

    template<class PlainScalar, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>
    Diff<Vector<PlainScalar>, DiffMode::Reverse, Order>::random_any(size_t len, Distribution& dist, RandomGenerator& gen) {
        return This(VectorType::random_any(len, dist, gen));
    }
}

namespace Physica {
    using namespace Core;

    template<class T, unsigned int Order>
    class Traits<Diff<Vector<T>, DiffMode::Reverse, Order>> : public Traits<Vector<T>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, DiffMode::Reverse, Order>;
    };
}
