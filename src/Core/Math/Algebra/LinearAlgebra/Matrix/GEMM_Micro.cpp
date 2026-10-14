/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixProduct/GEMM_Micro.h"

using namespace Physica;

namespace Physica::Internal {
    template<Scalar T, int LhsMajor, int RhsMajor>
    void gemm2x2(const T* lhs, const T* rhs, T* target) noexcept {
        using Pack = BestPacket<T, 4>::Type;
        static_assert(Pack::size() == 4, "[Error]: assign_micro() needs a 4 element packet");

        Pack a{}, b{};
        a.load(lhs);
        b.load(rhs);
        if constexpr (LhsMajor == MatrixMajor::Row)
            a = a.template permute<0, 2, 1, 3>();
        if constexpr (RhsMajor == MatrixMajor::Row)
            b = b.template permute<0, 2, 1, 3>();
        Pack c = fma(a.template permute<0, 1, 0, 1>(),
                     b.template permute<0, 0, 2, 2>(),
                     a.template permute<2, 3, 2, 3>() * b.template permute<1, 1, 3, 3>());
        c.store(target);
    }

    template<Scalar T, int LhsMajor, int RhsMajor>
    void gemm2x2_add(const T* lhs, const T* rhs, T* target) noexcept {
        using Pack = BestPacket<T, 4>::Type;
        static_assert(Pack::size() == 4, "[Error]: assign_micro() needs a 4 element packet");

        Pack a{}, b{}, c{};
        a.load(lhs);
        b.load(rhs);
        c.load(target);
        if constexpr (LhsMajor == MatrixMajor::Row)
            a = a.template permute<0, 2, 1, 3>();
        if constexpr (RhsMajor == MatrixMajor::Row)
            b = b.template permute<0, 2, 1, 3>();
        c = fma(a.template permute<0, 1, 0, 1>(), b.template permute<0, 0, 2, 2>(), c);
        c = fma(a.template permute<2, 3, 2, 3>(), b.template permute<1, 1, 3, 3>(), c);
        c.store(target);
    }

    template void gemm2x2<float64, MatrixMajor::Col, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    template void gemm2x2<float64, MatrixMajor::Row, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    template void gemm2x2<float64, MatrixMajor::Col, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;
    template void gemm2x2<float64, MatrixMajor::Row, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;

    template void gemm2x2<float32, MatrixMajor::Col, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    template void gemm2x2<float32, MatrixMajor::Row, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    template void gemm2x2<float32, MatrixMajor::Col, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;
    template void gemm2x2<float32, MatrixMajor::Row, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;

    template void gemm2x2_add<float64, MatrixMajor::Col, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    template void gemm2x2_add<float64, MatrixMajor::Row, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    template void gemm2x2_add<float64, MatrixMajor::Col, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;
    template void gemm2x2_add<float64, MatrixMajor::Row, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;

    template void gemm2x2_add<float32, MatrixMajor::Col, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    template void gemm2x2_add<float32, MatrixMajor::Row, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    template void gemm2x2_add<float32, MatrixMajor::Col, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;
    template void gemm2x2_add<float32, MatrixMajor::Row, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;
}
