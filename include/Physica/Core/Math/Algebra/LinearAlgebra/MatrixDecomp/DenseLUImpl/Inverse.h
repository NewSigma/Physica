/*
 * Copyright 2025-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    class Inverse<DenseLU<T, Pivot>> : public RValueMatrix<Inverse<DenseLU<T, Pivot>>> {
        using LU = DenseLU<T, Pivot>;
        using This = Inverse<LU>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::Tm;
    private:
        const LU& lu;
    public:
        Inverse(const LU& lu);
        Inverse(const This&) = default;
        Inverse(This&&) noexcept = default;
        ~Inverse() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t, size_t) const { noImpl(__func__); }

        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const LU& getDenseLU() const noexcept { return lu; }
        [[nodiscard]] size_t getRow() const noexcept { return lu.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return getRow(); }
    };

    template<Scalar T, bool Pivot>
    Inverse<DenseLU<T, Pivot>>::Inverse(const LU& lu) : lu(lu) {}

    template<Scalar T, bool Pivot>
    void Inverse<DenseLU<T, Pivot>>::assign(Matrix auto& target) const {
        target.assert_assign(*this);
        if constexpr (Internal::EnableMKL<typename LU::WorkingMatrix, decltype(target)>::value && Pivot)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Scalar T, bool Pivot>
    void Inverse<DenseLU<T, Pivot>>::assign_base(Matrix auto& target) const {
        static_assert(!Pivot, "[Error]: No impl");
        const auto& matrixLU = lu.getMatrixLU();
        matrixLU.tril_unit().inv().assign(target);
        (matrixLU.triu().inv() * target).assign(target);
    }
}

namespace Physica {
    template<Scalar T, bool P>
    class Traits<Inverse<DenseLU<T, P>>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::Col;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;

        using ExprType = DenseLU<T, P>;
    };
}

#ifdef PHYSICA_MKL
    #include "Inverse_MKL.h"
#endif
#include "InverseGEMV.h"
