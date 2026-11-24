/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        [[nodiscard]] T calc(size_t, size_t) const { noImpl(); }

        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const LU& getDenseLU() const noexcept { return lu; }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return lu.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Scalar T, bool Pivot>
    Inverse<DenseLU<T, Pivot>>::Inverse(const LU& lu) : lu(lu) {}

    template<Scalar T, bool Pivot>
    void Inverse<DenseLU<T, Pivot>>::assign(Matrix auto& target) const {
        target.assert_assign(*this);
        constexpr bool UseMKL = HasMKL() && Internal::EnableLAPACK<typename LU::WorkingMatrix, decltype(target)>::value && target.getMinorStrideAtCompile() == 1;
        if constexpr (UseMKL && Pivot)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Scalar T, bool Pivot>
    void Inverse<DenseLU<T, Pivot>>::assign_base(Matrix auto& target) const {
        const auto& matrixLU = lu.getMatrixLU();
        matrixLU.tril_unit().inv().assign(target);
        if constexpr (Pivot)
            target *= lu.getPerm().inv();
        (matrixLU.triu().inv() * target).assign(target);
    }
}

namespace Physica {
    template<Scalar T, bool P>
    class Traits<Inverse<DenseLU<T, P>>> {
    public:
        using ScalarType = T;
        using ExprType = DenseLU<T, P>;
    };
}

#ifdef PHYSICA_MKL
    #include "Inverse_MKL.h"
#endif
#include "InverseGEMV.h"
