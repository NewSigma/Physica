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
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"
#include "Physica/Core/Exception/NotImplementedException.h"

namespace Physica::Core {
    template<class ScalarType, bool IsSmallCell> class ChebyshevRadial;

    namespace Internal {
        template<class T, bool IsSmallCell>
        class Traits<ChebyshevRadial<T, IsSmallCell>> {
        public:
            using ScalarType = T;
        };
    }
    /**
     * Chebyshev descriptor introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 157, 114801 (2022); https://doi.org/10.1063/5.0106617
     */
    template<class ScalarType, bool IsSmallCell>
    class ChebyshevRadial : protected PairModel<ChebyshevRadial<ScalarType, IsSmallCell>> {
        using Base = PairModel<ChebyshevRadial<ScalarType, IsSmallCell>>;
        using PlainScalar = typename ScalarType::PlainScalar;
    public:
        using MDCellType = MDCell<PlainScalar, 3>;
        using ParticleType = typename MDCellType::ParticleType;
        using MassTypeMap = typename MDCellType::MassTypeMap;
        using DescriptorMatrix = DenseMatrix<ScalarType>;
        using DescriptorArray = Utils::Array<DescriptorMatrix>;
    private:
        MassTypeMap massTypeMap;
        unsigned int maxOrder;
    public:
        ChebyshevRadial(MassTypeMap massTypeMap_, unsigned int maxOrder_, ScalarType cutoff);
        ChebyshevRadial(const ChebyshevRadial&) = default;
        ChebyshevRadial(ChebyshevRadial&&) noexcept = default;
        ~ChebyshevRadial() = default;
        /* Operators */
        ChebyshevRadial& operator=(ChebyshevRadial obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] DescriptorArray project(const MDCellType& cell) const;
        void swap(ChebyshevRadial& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] unsigned int getMaxOrder() const noexcept { return maxOrder; }
    private:
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        /* Static members */
        [[nodiscard]] inline static ScalarType cutoff_functor(ScalarType r) noexcept;
        /* Friends */
        friend class PairModel<ChebyshevRadial<ScalarType, IsSmallCell>>;
    };

    template<class ScalarType, bool IsSmallCell>
    ChebyshevRadial<ScalarType, IsSmallCell>::ChebyshevRadial(MassTypeMap massTypeMap_, unsigned int maxOrder_, ScalarType cutoff)
            : Base(cutoff)
            , massTypeMap(std::move(massTypeMap_))
            , maxOrder(maxOrder_) {
        assert(maxOrder > 1 && "[Error]: Max order is too small to describe meaningful physics");
    }

    template<class ScalarType, bool IsSmallCell>
    void ChebyshevRadial<ScalarType, IsSmallCell>::swap(ChebyshevRadial& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        massTypeMap.swap(obj.massTypeMap);
        std::swap(maxOrder, obj.maxOrder);
    }

    template<class ScalarType, bool IsSmallCell>
    typename ChebyshevRadial<ScalarType, IsSmallCell>::DescriptorArray ChebyshevRadial<ScalarType, IsSmallCell>::project(
            const MDCellType& cell) const {
        const size_t numType = massTypeMap.size();
        DescriptorArray result(cell.getNumParticle(), maxOrder, numType, ScalarType(0));
        auto kernel = [this, maxOrder, &cell, &result](size_t i, size_t j, Vector3D r, ScalarType norm1, ScalarType norm2) {
            const ParticleType typeI = massTypeMap[cell.getMass(i)];
            const ParticleType typeJ = massTypeMap[cell.getMass(j)];
            auto descriptorDi = result[i].col(typeJ);
            auto descriptorDj = result[j].col(typeI);

            const ScalarType normalR = norm1 / Base::getCutoff();
            const ScalarType x0 = square(normalR - ScalarType(1));
            const ScalarType cutoffFactor = cutoff_functor(normalR);
            /* Make order 0 and 1 */ {
                descriptorDi[0] += cutoffFactor;
                descriptorDj[0] += cutoffFactor;

                const ScalarType descriptorD = x0 * cutoffFactor;
                descriptorDi[1] += descriptorD;
                descriptorDj[1] += descriptorD;
            }
            const ScalarType x = x0 * ScalarType(2) - ScalarType(1);
            ScalarType chebyshevN_2 = ScalarType(1);
            ScalarType chebyshevN_1 = x;
            ScalarType chebyshevN;
            for (size_t order = 2; order < maxOrder; ++order) {
                chebyshevN = ScalarType(2) * x * chebyshevN_1 - chebyshevN_2;
                const ScalarType descriptorD = chebyshevN * cutoffFactor + cutoffFactor;
                descriptorDi[order] += descriptorD;
                descriptorDj[order] += descriptorD;
                chebyshevN_2.swap(chebyshevN_1);
                chebyshevN_1.swap(chebyshevN);
            }
        };
        Base::forPairInCutoff(cell.getLattice(), cell.getPos(), kernel);
        return result;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType ChebyshevRadial<ScalarType, IsSmallCell>::pot_functor(size_t, size_t, ScalarType, ScalarType) const {
        throw NotImplementedException("[Error]: This function is disabled");
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType ChebyshevRadial<ScalarType, IsSmallCell>::force_functor(size_t, size_t, ScalarType, ScalarType) const {
        throw NotImplementedException("[Error]: This function is disabled");
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType ChebyshevRadial<ScalarType, IsSmallCell>::cutoff_functor(ScalarType normalR) noexcept {
        assert(normalR.isPositive() && (normalR <= ScalarType(1)) && "[Error]: Distance out of cutoff");
        return square(cos(normalR * ScalarType(M_PI_2)));
    }
}
