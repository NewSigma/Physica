/*
 * Copyright 2024 WeiBo He.
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

#include "Physica/Core/MultiPrecision/Differentiable.h"
#include "Vector.h"

namespace Physica::Core {
    namespace Internal {
        template<class T>
        class Traits<Differentiable<Vector<T>, DiffMode::Reverse>> : public Traits<Vector<T>> {
        public:
            using ScalarType = Differentiable<T, DiffMode::Reverse>;
        };
    }

    template<class PlainScalar>
    class Differentiable<Vector<PlainScalar>, DiffMode::Reverse>
            : public RValueVector<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");

        using VectorType = Vector<PlainScalar>;
        using This = Differentiable<VectorType, DiffMode::Reverse>;
        using Base = RValueVector<This>;
        using TracerType = DiffTracer<PlainScalar>;
        using SegmentType = TraceSegment<PlainScalar>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Differentiable();
        Differentiable(size_t length);
        Differentiable(VectorType values);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(const Differentiable&) = default;
        Differentiable& operator=(Differentiable&&) noexcept = default;
        /* Operations */
        template<class RandomGenerator> inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator> inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(Differentiable& obj) noexcept { std::swap(*this, obj); }
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

    template<class PlainScalar>
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::Differentiable()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<class PlainScalar>
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::Differentiable(size_t length)
            : traceSeg(TracerType::getInstance().pushSegment(length)) {}

    template<class PlainScalar>
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::Differentiable(VectorType values)
            : traceSeg(TracerType::getInstance().pushSegment(std::move(values))) {}

    template<class PlainScalar>
    template<class RandomGenerator>
    inline void Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getLength(), gen);
    }

    template<class PlainScalar>
    template<class RandomGenerator>
    inline void Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getLength(), gen);
    }

    template<class PlainScalar>
    template<class Distribution, class RandomGenerator>
    inline void Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getLength(), dist, gen);
    }

    template<class PlainScalar>
    inline typename Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::ScalarType
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        return traceSeg[index];
    }

    template<class PlainScalar>
    template<class RandomGenerator>
    inline Differentiable<Vector<PlainScalar>, DiffMode::Reverse>
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::random_uniform(size_t len, RandomGenerator& gen) {
        return This(VectorType::random_uniform(len, gen));
    }

    template<class PlainScalar>
    template<class RandomGenerator>
    inline Differentiable<Vector<PlainScalar>, DiffMode::Reverse>
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::random_normal(size_t len, RandomGenerator& gen) {
        return This(VectorType::random_normal(len, gen));
    }

    template<class PlainScalar>
    template<class Distribution, class RandomGenerator>
    inline Differentiable<Vector<PlainScalar>, DiffMode::Reverse>
    Differentiable<Vector<PlainScalar>, DiffMode::Reverse>::random_any(size_t len, Distribution& dist, RandomGenerator& gen) {
        return This(VectorType::random_any(len, dist, gen));
    }
}
