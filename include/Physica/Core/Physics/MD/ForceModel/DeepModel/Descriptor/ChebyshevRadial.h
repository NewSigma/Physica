/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"
#include "Physica/Core/Exception/NoImplException.h"

namespace Physica::Core {
    /**
     * Chebyshev descriptor introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 157, 114801 (2022); https://doi.org/10.1063/5.0106617
     */
    template<Scalar T, bool IsSmallCell>
    class ChebyshevRadial : protected PairModel<ChebyshevRadial<T, IsSmallCell>> {
        using Base = PairModel<ChebyshevRadial<T, IsSmallCell>>;
        using ValueType = T::ValueType;
    public:
        using MDCellType = MDCell<ValueType, 3>;
        using ParticleType = MDCellType::ParticleType;
        using MassTypeMap = MDCellType::MassTypeMap;
        using DescriptorMatrix = DenseMatrix<T>;
        using DescriptorArray = Array<DescriptorMatrix>;
    private:
        MassTypeMap massTypeMap;
        unsigned int maxOrder;
    public:
        ChebyshevRadial(MassTypeMap massTypeMap_, unsigned int maxOrder_, T cutoff);
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
        [[nodiscard]] inline T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] inline T force_functor(size_t i, size_t j, T r, T r2) const;
        /* Static members */
        [[nodiscard]] inline static T cutoff_functor(T r) noexcept;
        /* Friends */
        friend class PairModel<ChebyshevRadial<T, IsSmallCell>>;
    };

    template<Scalar T, bool IsSmallCell>
    ChebyshevRadial<T, IsSmallCell>::ChebyshevRadial(MassTypeMap massTypeMap_, unsigned int maxOrder_, T cutoff)
            : Base(cutoff)
            , massTypeMap(std::move(massTypeMap_))
            , maxOrder(maxOrder_) {
        assert(maxOrder > 1 && "[Error]: Max order is too small to describe meaningful physics");
    }

    template<Scalar T, bool IsSmallCell>
    void ChebyshevRadial<T, IsSmallCell>::swap(ChebyshevRadial& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        massTypeMap.swap(obj.massTypeMap);
        std::swap(maxOrder, obj.maxOrder);
    }

    template<Scalar T, bool IsSmallCell>
    ChebyshevRadial<T, IsSmallCell>::DescriptorArray ChebyshevRadial<T, IsSmallCell>::project(
            const MDCellType& cell) const {
        const size_t numType = massTypeMap.size();
        DescriptorArray result(cell.getNumParticle(), maxOrder, numType, T(0));
        auto kernel = [this, &cell, &result](size_t i, size_t j, Vector3D<T> r, T norm1, T norm2) {
            const ParticleType typeI = massTypeMap[cell.getMass(i)];
            const ParticleType typeJ = massTypeMap[cell.getMass(j)];
            auto descriptorDi = result[i].col(typeJ);
            auto descriptorDj = result[j].col(typeI);

            const T normalR = norm1 / Base::getCutoff();
            const T x0 = square(normalR - T(1));
            const T cutoffFactor = cutoff_functor(normalR);
            /* Make order 0 and 1 */ {
                descriptorDi[0] += cutoffFactor;
                descriptorDj[0] += cutoffFactor;

                const T descriptorD = x0 * cutoffFactor;
                descriptorDi[1] += descriptorD;
                descriptorDj[1] += descriptorD;
            }
            const T x = x0 * T(2) - T(1);
            T chebyshevN_2 = T(1);
            T chebyshevN_1 = x;
            T chebyshevN;
            for (size_t order = 2; order < maxOrder; ++order) {
                chebyshevN = T(2) * x * chebyshevN_1 - chebyshevN_2;
                const T descriptorD = chebyshevN * cutoffFactor + cutoffFactor;
                descriptorDi[order] += descriptorD;
                descriptorDj[order] += descriptorD;
                chebyshevN_2.swap(chebyshevN_1);
                chebyshevN_1.swap(chebyshevN);
            }
        };
        Base::forPairInCutoff(cell.getLattice(), cell.getPos(), kernel);
        return result;
    }

    template<Scalar T, bool IsSmallCell>
    inline T ChebyshevRadial<T, IsSmallCell>::pot_functor(size_t, size_t, T, T) const {
        noImpl("[Error]: This function is disabled");
    }

    template<Scalar T, bool IsSmallCell>
    inline T ChebyshevRadial<T, IsSmallCell>::force_functor(size_t, size_t, T, T) const {
        noImpl("[Error]: This function is disabled");
    }

    template<Scalar T, bool IsSmallCell>
    inline T ChebyshevRadial<T, IsSmallCell>::cutoff_functor(T normalR) noexcept {
        assert(normalR.isPositive() && (normalR <= T(1)) && "[Error]: Distance out of cutoff");
        return square(cos(normalR * T(M_PI_2)));
    }
}

namespace Physica {
    template<class T, bool IsSmallCell>
    class Traits<Core::ChebyshevRadial<T, IsSmallCell>> {
    public:
        using ScalarType = T;
    };
}
