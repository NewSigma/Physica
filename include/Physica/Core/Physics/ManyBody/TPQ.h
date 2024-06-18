/*
 * Copyright 2024 WeiBo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h>
#include "LatticeHamilton.h"

namespace Physica::Core {
    /**
     * NVE ensemble referenced from [1]
     * NVT ensemble referenced from [2]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 108, 240401 (2012); https://doi.org/10.1103/PhysRevLett.108.240401
     * [2] Phys. Rev. Lett. 111, 010401 (2013); https://doi.org/10.1103/PhysRevLett.108.240401
     */
    template<class ScalarType>
    class TPQ : public Vector<ScalarType> {
        using This = TPQ<ScalarType>;
        using Base = Vector<ScalarType>;
        using RealType = typename ScalarType::RealType;

        RealType beta;
    public:
        TPQ(const This&) = default;
        TPQ(This&&) noexcept = default;
        ~TPQ() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class ModelType>
        void nvt_step(const LatticeHamilton<ModelType>& hamiltonH, RealType deltaBeta);
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] static This random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static This random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] static This random_any(size_t len, Distribution& dist, RandomGenerator& gen);
    private:
        TPQ(Base v_);
    };

    template<class ScalarType>
    TPQ<ScalarType>::TPQ(Base v_) : Base(std::move(v_)), beta(0) {
        Base::toUnit();
    }

    template<class ScalarType>
    template<class ModelType>
    void TPQ<ScalarType>::nvt_step(const LatticeHamilton<ModelType>& hamiltonH, RealType deltaBeta) {
        Vector<ScalarType> dot = exp((-deltaBeta * 0.5) * hamiltonH.getDerived()) * (*this);
        Base::swap(dot);
        beta += deltaBeta;
    }

    template<class ScalarType>
    void TPQ<ScalarType>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(*this);
        beta.swap(obj.beta);
    }

    template<class ScalarType>
    template<class RandomGenerator>
    TPQ<ScalarType> TPQ<ScalarType>::random_uniform(size_t len, RandomGenerator& gen) {
        return This(Base::random_uniform(len, gen));
    }

    template<class ScalarType>
    template<class RandomGenerator>
    TPQ<ScalarType> TPQ<ScalarType>::random_normal(size_t len, RandomGenerator& gen) {
        return This(Base::random_normal(len, gen));
    }

    template<class ScalarType>
    template<class Distribution, class RandomGenerator>
    TPQ<ScalarType> TPQ<ScalarType>::random_any(size_t len, Distribution& dist, RandomGenerator& gen) {
        return This(Base::random_any(len, dist, gen));
    }
}
