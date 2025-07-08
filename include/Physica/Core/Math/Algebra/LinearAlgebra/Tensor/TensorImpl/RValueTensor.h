/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/CRTPBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"
#include "TensorBase.h"

namespace Physica {
    template<class Derived> class LValueTensor;
    template<class> class FlattenL;
    template<class T> class RealTensor;
    template<class T> class ImagTensor;
    template<class T> class SquaredNormTensor;
    template<class T> class NormTensor;
    template<class T> class ValueTensor;
    template<class T, int GradOrder> class GradTensor;

    template<class Derived>
    class RValueTensor : public CRTPBase<RValueTensor<Derived>>, public TensorBase {
        using This = RValueTensor<Derived>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
        constexpr static int Dim = Traits<Derived>::Dim;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
    protected:
        using T = ScalarType;
        using IndexArray = Array<size_t, Dim>;
    public:
        /* Operations */
        void assign(Tensor auto& x) const;

        [[nodiscard]] auto calc(size_t x, size_t y, size_t z) const { return calc({x, y, z}); }
        [[nodiscard]] auto calc(Index3D index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] size_t toIndex1D(const IndexArray& indices) const noexcept;
        [[nodiscard]] IndexArray toIndexND(size_t index) const noexcept;

        void forND(std::invocable<T, IndexArray> auto fn) const;

        auto reals() const noexcept;
        auto imags() const noexcept;
        auto squaredNorms() const noexcept;
        auto norms() const noexcept;
        auto values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] auto getShape(int dim) const noexcept { return Base::getDerived().getShape(dim); }
        [[nodiscard]] IndexArray getShape() const noexcept;
        [[nodiscard]] auto getDim() const noexcept;
        [[nodiscard]] size_t getSize() const noexcept;
        /* Static members */
        using TensorBase::forPointIndexInTensor;
        template<Scalar T, bool IsUnitLattice>
        static void forPointIndexInTensor(const RValueTensor& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, std::invocable<Vector3D<T>, Index3D> auto fn) {
            forPointIndexInTensor<T, IsUnitLattice>(grid.getDim(), lattice, fn); // FIXME: NVCC 12.8 rejects valid if we put it in impl file
        }

        [[nodiscard]] static size_t toSize(const IndexArray& shape);
    protected:
        template<int GradOrder>
        auto grads_impl() const noexcept;
    };
}

namespace Physica {
    template<class T>
    class Traits<RValueTensor<T>> {
    public:
        using Derived = T;
    };
}

#include "RValueTensorImpl/RValueTensorImpl.h"
#include "RValueTensorImpl/TensorConvert.h"
#include "TensorExpr.h"
