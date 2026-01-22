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
        [[nodiscard]] T calc(size_t, size_t) const { noImpl(__func__); }

        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const LU& getDenseLU() const noexcept { return lu; }
        [[nodiscard]] size_t getRow() const noexcept { return lu.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return getRow(); }
    };

    template<Scalar T, bool Pivot>
    device_obj<Inverse<DenseLU<T, Pivot>>>::device_obj(const LU& lu) : lu(lu) {}

    template<Scalar T, bool Pivot>
    void device_obj<Inverse<DenseLU<T, Pivot>>>::assign(Matrix auto& target) const {
        static_assert(!Pivot, "[Error]: No impl");
        target.assert_assign(*this);
        
        const auto& matrixLU = lu.getMatrixLU();
        matrixLU.tril_unit().inv().assign(target);
        (matrixLU.triu().inv() * target).assign(target);
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
