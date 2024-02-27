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
#include "DenseMatrix.h"

namespace Physica::Core {
    namespace Internal {
        template<class T, int Option, unsigned int Order>
        class Traits<Differentiable<DenseMatrix<T, Option>, DiffMode::Reverse, Order>> : public Traits<DenseMatrix<T, Option>> {
        public:
            using ScalarType = Differentiable<T, DiffMode::Reverse, Order>;
        };
    }

    template<class PlainScalar, int Option, unsigned int Order>
    class Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
            : public RValueMatrix<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
            , public DenseMatrixDim<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>, Dynamic, Dynamic, Dynamic, Dynamic> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using PlainMatrix = DenseMatrix<PlainScalar, Option>;
        using This = Differentiable<PlainMatrix, DiffMode::Reverse, Order>;
        using Base = RValueMatrix<This>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic, Dynamic, Dynamic>;
        using TracerType = DiffTracer<PlainScalar, Order>;
        using SegmentType = typename TracerType::SegmentType;
    public:
        using device_obj_type = device_obj<This>;
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Differentiable();
        Differentiable(size_t row, size_t column);
        Differentiable(PlainMatrix values);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(const Differentiable&) = default;
        Differentiable& operator=(Differentiable&&) noexcept = default;
        /* Operations */
        [[nodiscard]] inline ScalarType calc(size_t row, size_t col) const;

        template<class RandomGenerator> inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator> inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(Differentiable& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getRow;
        using Dim::getColumn;
        [[nodiscard]] size_t getSize() const noexcept { return traceSeg.getLength(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t column, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t row, size_t column, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t row, size_t column, Distribution& dist, RandomGenerator& gen);
    private:
        friend class device_obj<This>;
    };

    template<class PlainScalar, int Option, unsigned int Order>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::Differentiable()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<class PlainScalar, int Option, unsigned int Order>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::Differentiable(size_t row, size_t column)
            : traceSeg(TracerType::getInstance().pushSegment(row * column)) {}

    template<class PlainScalar, int Option, unsigned int Order>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::Differentiable(PlainMatrix values)
            : traceSeg(TracerType::getInstance().pushSegment(values.flatten())) {}

    template<class PlainScalar, int Option, unsigned int Order>
    inline typename Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::ScalarType
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::calc(size_t row, size_t col) const {
        if constexpr (Base::isRowMatrix)
            return traceSeg[row * getColumn() + col];
        else
            return traceSeg[col * getRow() + row];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline void Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getRow(), getColumn(), gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline void Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getRow(), getColumn(), gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline void Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getRow(), getColumn(), dist, gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_uniform(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(PlainMatrix::random_uniform(row, column, gen));
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_normal(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(PlainMatrix::random_normal(row, column, gen));
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_any(
            size_t row, size_t column, Distribution& dist, RandomGenerator& gen) {
        return This(PlainMatrix::random_any(row, column, dist, gen));
    }
}
