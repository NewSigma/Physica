/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"

namespace Physica {
    template<Scalar> class LUMatrixL;

    template<Scalar T, bool Pivot>
    class LUDecomp {
        using This = LUDecomp;
        using WorkingMatrix = DenseMatrix<T>;
        using BiasArray = std::conditional<Pivot, PermMatrix<T>, PlainStruct<void>>::type;

        constexpr static bool isComplex = T::isComplex;
        using Tc = T::ComplexType;
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
    private:
        WorkingMatrix working;
        [[no_unique_address]] BiasArray perm;
    public:
        LUDecomp() = default;
        LUDecomp(size_t size);
        LUDecomp(const Matrix auto& source);
        LUDecomp(const This&) = default;
        LUDecomp(This&&) noexcept = default;
        ~LUDecomp() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void compute_mkl(const Matrix auto& source);
        void compute_base(const Matrix auto& source);

        [[nodiscard]] T det() const noexcept;

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] const auto& getMatrixLU() const noexcept { return working; }
        [[nodiscard]] auto getMatrixL() const noexcept { return LUMatrixL<T>(*this); }
        [[nodiscard]] auto getMatrixU() const noexcept { return working.triu(); }
        [[nodiscard]] const PermMatrix<T>& getPerm() const noexcept;
    private:
        void pre_compute(const Matrix auto& source) const noexcept;
        void decomp_col(size_t col);
    };
}

#ifdef PHYSICA_MKL
    #include "LUDecompImpl/LUDecompImpl.h"
    #include "LUDecompImpl/LUDecomp_MKL.h"
#endif
