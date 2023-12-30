/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"

namespace Physica::Core {
    template<class ScalarType, bool IsSmallCell> class GaussKernel;

    namespace Internal {
        template<class T, bool IsSmallCell>
        class Traits<GaussKernel<T, IsSmallCell>> {
        public:
            using ScalarType = T;
        };
    }
    /**
     * Reference:
     * [1] Phys. Rev. Lett. 98, 146401 (2007); https://doi.org/10.1103/PhysRevLett.98.146401
     */
    template<class ScalarType, bool IsSmallCell>
    class GaussKernel : protected PairModel<GaussKernel<ScalarType, IsSmallCell>> {
        using Base = PairModel<GaussKernel<ScalarType, IsSmallCell>>;
        using VectorType = Vector<ScalarType>;
        using CellType = PeriodicCell<ScalarType, 3>; // 3 dim because low dimension is seldom used.

        ScalarType paramEta;
        ScalarType distR;
    public:
        GaussKernel(ScalarType paramEta_, ScalarType distR_, ScalarType cutoff);
        GaussKernel(const GaussKernel&) = default;
        GaussKernel(GaussKernel&&) noexcept = default;
        ~GaussKernel() = default;
        /* Operators */
        GaussKernel& operator=(GaussKernel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorType calc(const CellType& cell) const;
        void swap(GaussKernel& __restrict obj) noexcept;
    private:
        using Base::force_functor;
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType cutoff_functor(ScalarType r) const;
        /* Friends */
        friend class PairModel<GaussKernel<ScalarType, IsSmallCell>>;
    };

    template<class ScalarType, bool IsSmallCell>
    GaussKernel<ScalarType, IsSmallCell>::GaussKernel(ScalarType paramEta_, ScalarType distR_, ScalarType cutoff)
        : Base(cutoff)
        , paramEta(std::move(paramEta_))
        , distR(std::move(distR_)) {}

    template<class ScalarType, bool IsSmallCell>
    void GaussKernel<ScalarType, IsSmallCell>::swap(GaussKernel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        paramEta.swap(obj.paramEta);
        distR.swap(obj.distR);
    }

    template<class ScalarType, bool IsSmallCell>
    typename GaussKernel<ScalarType, IsSmallCell>::VectorType GaussKernel<ScalarType, IsSmallCell>::calc(const CellType& cell) const {
        VectorType result(cell.getNumParticle(), 0);
        auto kernel = [this, &result](size_t i, size_t j, Vector3D r, ScalarType norm1, ScalarType norm2) {
            const ScalarType value = exp(-paramEta * square(norm1 - distR)) * cutoff_functor(norm1);
            result[i] += value;
            result[j] += value;
        };
        Base::forPairInCutoff(cell.getLattice(), cell.getPos(), kernel);
        return result;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType GaussKernel<ScalarType, IsSmallCell>::pot_functor(
            [[maybe_unused]] size_t i,
            [[maybe_unused]] size_t j,
            [[maybe_unused]] ScalarType r,
            [[maybe_unused]] ScalarType r2) const {
        return ScalarType(0);
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType GaussKernel<ScalarType, IsSmallCell>::cutoff_functor(ScalarType r) const {
        if (r > Base::getCutoff())
            return ScalarType(0);
        return (cos(ScalarType(M_PI) * r / Base::getCutoff()) + ScalarType(1)) * 0.5;
    }
}
