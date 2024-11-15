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

#include "PairModel.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Physics/MD/MDCell.cuh"
#include "Physica/Core/Physics/MD/MDImpl/CellList.cuh"
#include "Physica/Core/Utils/Allocator/PageLockedAllocator.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<PairModel<Derived>> : public CRTPBase<device_obj<PairModel<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: Nested device_obj is unnecessary");
        using host_obj = PairModel<Derived>;
        using This = device_obj<PairModel<Derived>>;
        using Base = CRTPBase<This>;
        using TraitType = Traits<Derived>;

        constexpr static bool IsPotDependOnAtomIndex = TraitType::IsPotDependOnAtomIndex;
        constexpr static bool IsSmallCell = TraitType::IsSmallCell;
    public:
        using ScalarType = typename TraitType::ScalarType;
        constexpr static int Dim = host_obj::Dim;
        constexpr static int NumVirialElem = Dim * Dim;

        using ValueType = typename ScalarType::ValueType;
        using MDCellType = MDCell<ScalarType>;
        using DeviceMDCell = device_obj<MDCellType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using InvLatticeMatrix = typename MDCellType::InvLatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using CellListType = CellList<ScalarType>;
        using DeviceCellList = device_obj<CellListType>;
        using Index3D = typename GridBase::Index3D;
        using Vector3D = device_obj<Vector3D<ScalarType>>;
        using ForceBufferType = device_obj<DenseMatrix<ScalarType>>;
        using VirialBufferType = device_obj<DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Element, NumVirialElem>>;
        using PageLockedVector = DenseVector<ScalarType, Dynamic, PageLockedAllocator<ScalarType>>;
    private:
        ScalarType cutoff;
        ScalarType squared_cutoff;
        ScalarType pot_shift;
        DeviceMDCell cell;
        DeviceCellList cellList;
        ForceBufferType forceBuffer;
        VirialBufferType virialBuffer;
        PageLockedVector swapBuffer;
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] __host__ __device__ inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] __host__ __device__ inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;

        [[nodiscard]] ScalarType potentialV(const MDCellType& hostCell) const;

        template<class Executor>
        [[nodiscard]] VectorND<ScalarType> force(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos);
        template<class Executor>
        [[nodiscard]] inline VectorND<ScalarType> force(const MDCellType& hostCell);

        template<class VectorType, class Executor>
        void forceAsync(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos,
                ContinuousVector<VectorType>& result);
        template<class VectorType, class Executor>
        inline void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result);
        template<class Executor>
        [[nodiscard]] inline VectorND<ScalarType> force_short(const MDCellType& cell);
        template<class Executor>
        [[nodiscard]] inline VectorND<ScalarType> force_long(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos);
        [[nodiscard]] inline LatticeMatrix virial(const MDCellType& hostCell);
        void swap(device_obj& __restrict obj) noexcept;

        __device__ void forceKernelImpl();
        __device__ void postForceKernelImpl();
        __device__ void virialKernelImpl();
        __device__ void postVirialKernelImpl();
        /* Getters */
        [[nodiscard]] __host__ __device__ const ScalarType& getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] __host__ __device__ const ScalarType& getSquaredCutoff() const noexcept { return squared_cutoff; }
        [[nodiscard]] __device__ const DeviceMDCell& getCell() const noexcept { return cell; }
        /* Setters */
        void setCutoff(ScalarType cutoff_);
    protected:
        device_obj() = default;
        device_obj(size_t numParticle);
        device_obj(size_t numParticle, ScalarType cutoff_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void preParallel(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos,
                dim3& gridDims,
                size_t& numThread);
        template<class Functor>
        [[nodiscard]] __device__ size_t forPairInCutoff(Functor func) const;
    };
}

namespace Physica {
    template<class T>
    class Traits<Core::device_obj<Core::PairModel<T>>> : public Traits<Core::PairModel<T>> {
    public:
        using Derived = Core::device_obj<T>;
    };
}

#include "PairModelImpl.cuh"
