/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"

namespace Physica {
    template<Scalar> class LUMatrixL;
    template<Scalar, bool Pivot> class DenseLUInv;

    template<Scalar T, bool Pivot = true>
    class DenseLU {
        using This = DenseLU;
        using Tr = T::RealType;
        using Tc = T::ComplexType;
        using PermType = std::conditional<Pivot, PermMatrix<Tr>, Empty>::type;
    public:
        using WorkingMatrix = DenseMatrix<T>;
    private:
        WorkingMatrix working;
        [[no_unique_address]] PermType perm;
    public:
        DenseLU() = default;
        DenseLU(size_t size);
        DenseLU(const Matrix auto& source);
        DenseLU(const This&) = default;
        DenseLU(This&&) noexcept = default;
        ~DenseLU() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute();
        void compute_mkl();
        void compute_base();

        void compute(const Matrix auto& source);
        void compute_mkl(const Matrix auto& source);
        void compute_base(const Matrix auto& source);

        [[nodiscard]] T det() const noexcept;
        [[nodiscard]] Tr lnAbsDet() const noexcept;
        [[nodiscard]] T sgndet() const noexcept;
        [[nodiscard]] auto inv() const noexcept;

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] constexpr size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return working.getOrder(); }
        [[nodiscard]] auto& getMatrixLU(this auto&& self) noexcept;
        [[nodiscard]] auto getMatrixL() const noexcept { return working.tril_unit(); }
        [[nodiscard]] auto getMatrixU() const noexcept { return working.triu(); }
        [[nodiscard]] const auto& getPerm() const noexcept;
        [[nodiscard]] consteval static bool isPivot() noexcept { return Pivot; }
    private:
        void decomp(Matrix auto& block, size_t offset);
    };
}

#include "DenseLUImpl/DenseLUImpl.h"
#ifdef PHYSICA_MKL
    #include "DenseLUImpl/DenseLU_MKL.h"
#endif
#include "DenseLUImpl/Inverse.h"
