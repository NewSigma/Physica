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
#include "DenseMatrix.h"

namespace Physica::Core {
    template<class PlainScalar, int Option, unsigned int Order>
    class Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
            : public RValueMatrix<Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
            , public DenseMatrixDim<Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>, Dynamic, Dynamic> {
        using PlainMatrix = DenseMatrix<PlainScalar, Option>;
        using This = Diff<PlainMatrix, DiffMode::Reverse, Order>;
        using Base = RValueMatrix<This>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic>;
        using TracerType = DiffTracer<PlainScalar, Order>;
        using SegmentType = typename TracerType::SegmentType;
    public:
        using device_obj_type = device_obj<This>;
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Diff();
        Diff(size_t row, size_t column);
        Diff(PlainMatrix values);
        Diff(const Diff&) = default;
        Diff(Diff&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        Diff& operator=(const Diff&) = default;
        Diff& operator=(Diff&&) noexcept = default;
        /* Operations */
        [[nodiscard]] inline ScalarType calc(size_t row, size_t col) const;

        template<class RandomGenerator>
        inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(Diff& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
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
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::Diff()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<class PlainScalar, int Option, unsigned int Order>
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::Diff(size_t row, size_t column)
            : traceSeg(TracerType::getInstance().pushSegment(row * column)) {}

    template<class PlainScalar, int Option, unsigned int Order>
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::Diff(PlainMatrix values)
            : traceSeg(TracerType::getInstance().pushSegment(values.flatten())) {}

    template<class PlainScalar, int Option, unsigned int Order>
    inline typename Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::ScalarType
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::calc(size_t row, size_t col) const {
        if constexpr (Base::isRowMatrix)
            return traceSeg[row * getColumn() + col];
        else
            return traceSeg[col * getRow() + row];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline void Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getRow(), getColumn(), gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline void Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getRow(), getColumn(), gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline void Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getRow(), getColumn(), dist, gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_uniform(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(PlainMatrix::random_uniform(row, column, gen));
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_normal(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(PlainMatrix::random_normal(row, column, gen));
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
    Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>::random_any(
            size_t row, size_t column, Distribution& dist, RandomGenerator& gen) {
        return This(PlainMatrix::random_any(row, column, dist, gen));
    }
}

namespace Physica {
    using namespace Core;

    template<class T, int Option, unsigned int Order>
    class Traits<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>> : public Traits<DenseMatrix<T, Option>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, DiffMode::Reverse, Order>;
    };
}
