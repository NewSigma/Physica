/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Scalar T, int Dim>
    class IntegrateRange {
        //We consider dimension which is larger than 3 is unsuitable to allocate domain data on stack.
        static_assert(Dim <= 3, "Dimension larger than 3 must be set Dynamic.");
    public:
        using VectorType = DenseVector<T, Dim>;
    private:
        VectorType range_from;
        VectorType range_to;
    public:
        IntegrateRange(VectorType range_from_, VectorType range_to_);
        /* Getters */
        [[nodiscard]] const VectorType& from() const { return range_from; }
        [[nodiscard]] const VectorType& to() const { return range_to; }
    };

    template<Scalar T, int Dim>
    IntegrateRange<T, Dim>::IntegrateRange(VectorType range_from_, VectorType range_to_)
            : range_from(std::move(range_from_)), range_to(std::move(range_to_)) {}
}
