/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    class PHYSICA_API Gnuplot {
    public:
        using DataArray = Array<VectorND<float64>>;
    private:
        DataArray xDatas;
        DataArray yDatas;
    public:
        Gnuplot() = default;
        Gnuplot(const Gnuplot&) = default;
        Gnuplot(Gnuplot&&) noexcept = default;
        ~Gnuplot() = default;
        /* Operators */
        Gnuplot& operator=(const Gnuplot&) = default;
        Gnuplot& operator=(Gnuplot&&) noexcept = default;
        friend std::ostream& operator<<(std::ostream& os, const Gnuplot& gnuplot);
        friend std::istream& operator>>(std::istream& is, Gnuplot& gnuplot);
        /* Getters */
        [[nodiscard]] const DataArray& getXDatas() const noexcept { return xDatas; }
        [[nodiscard]] const DataArray& getYDatas() const noexcept { return yDatas; }
    };

    PHYSICA_API std::ostream& operator<<(std::ostream& os, const Gnuplot& gnuplot);
    PHYSICA_API std::istream& operator>>(std::istream& is, Gnuplot& gnuplot);
}
