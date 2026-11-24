/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.cuh"

namespace Physica {
    template<Scalar T, bool Pivot>
    class device_obj<Inverse<DenseLU<T, Pivot>>> : public device_obj<RValueMatrix<Inverse<DenseLU<T, Pivot>>>> {
        using host_obj = Inverse<DenseLU<T, Pivot>>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using LU = device_obj<DenseLU<T, Pivot>>;
    public:
        using Base::isComplex;
    private:
        const LU& lu;
    public:
        device_obj(const LU& lu);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] T calc(size_t, size_t, instanceof_x<ThreadBlock> auto) const { noImpl(); }

        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const LU& getDenseLU() const noexcept { return lu; }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return lu.getOrder(); }
    };

    template<Scalar T, bool Pivot>
    device_obj<Inverse<DenseLU<T, Pivot>>>::device_obj(const LU& lu) : lu(lu) {}

    template<Scalar T, bool Pivot>
    void device_obj<Inverse<DenseLU<T, Pivot>>>::assign(Matrix auto& target) const {
        target.assert_assign(*this);
        if constexpr (Pivot) {
            target.zeros();
            target.diag() = T(1);
            lu.solve(target);
        }
        else {
            const auto& matrixLU = lu.getMatrixLU();
            matrixLU.tril_unit().inv().assign(target);
            (matrixLU.triu().inv() * target).assign(target);
        }
    }
}

namespace Physica {
    template<Scalar T, bool P>
    class Traits<device_obj<Inverse<DenseLU<T, P>>>> : public Traits<Inverse<DenseLU<T, P>>> {};
}

#ifdef PHYSICA_MKL
    #include "Inverse_MKL.h"
#endif
#include "InverseGEMV.h"
