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
        template<class T, int Option>
        class Traits<Differentiable<DenseMatrix<T, Option>, DiffMode::Reverse>> : public Traits<DenseMatrix<T, Option>> {
        public:
            using ScalarType = Differentiable<T, DiffMode::Reverse>;
        };
    }

    template<class PlainScalar, int Option>
    class Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>
            : public RValueMatrix<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>>
            , public DenseMatrixDim<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>, Dynamic, Dynamic, Dynamic, Dynamic> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using MatrixType = DenseMatrix<PlainScalar, Option>;
        using This = Differentiable<MatrixType, DiffMode::Reverse>;
        using Base = RValueMatrix<This>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic, Dynamic, Dynamic>;
        using TracerType = DiffTracer<PlainScalar>;
        using SegmentType = TraceSegment<PlainScalar>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Differentiable();
        Differentiable(MatrixType values);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(const Differentiable&) = default;
        Differentiable& operator=(Differentiable&&) noexcept = default;
        /* Operations */
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
    };

    template<class PlainScalar, int Option>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>::Differentiable()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<class PlainScalar, int Option>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>::Differentiable(MatrixType values)
            : traceSeg(TracerType::getInstance().pushSegment(values.flatten())) {}

    template<class PlainScalar, int Option>
    template<class RandomGenerator>
    inline Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>::random_uniform(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(MatrixType::random_uniform(row, column, gen));
    }

    template<class PlainScalar, int Option>
    template<class RandomGenerator>
    inline Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>::random_normal(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(MatrixType::random_normal(row, column, gen));
    }

    template<class PlainScalar, int Option>
    template<class Distribution, class RandomGenerator>
    inline Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>
    Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse>::random_any(
            size_t row, size_t column, Distribution& dist, RandomGenerator& gen) {
        return This(MatrixType::random_any(row, column, dist, gen));
    }
}
