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

#include "PairModel.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Physics/MD/MDCell.cuh"
#include "Physica/Core/Physics/MD/MDImpl/CellList.cuh"
#include "Physica/Core/Utils/Allocator/PageLockedAllocator.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<PairModel<Derived>> : public CRTPBase<device_obj<PairModel<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: Nested device_obj is unnecessary");
        using host_obj = PairModel<Derived>;
        using This = device_obj<PairModel<Derived>>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;

        constexpr static bool IsPotDependOnAtomIndex = TraitsType::IsPotDependOnAtomIndex;
        constexpr static bool IsSmallCell = TraitsType::IsSmallCell;
    public:
        constexpr static int Dim = host_obj::Dim;
        constexpr static int NumVirialElem = Dim * Dim;
        using ScalarType = TraitsType::ScalarType;
        using T = ScalarType;
        using Tv = T::ValueType;
        using MDCellType = MDCell<T>;
        using DeviceMDCell = device_obj<MDCellType>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using InvLatticeMatrix = MDCellType::InvLatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using CellListType = CellList<T>;
        using DeviceCellList = device_obj<CellListType>;
        using DeviceVector3D = device_obj<Vector3D<T>>;
        using ForceBufferType = device_obj<DenseMatrix<T>>;
        using VirialBufferType = device_obj<DenseMatrix<T, MatrixMajor::Col, NumVirialElem>>;
        using PageLockedVector = DenseVector<T, Dynamic, PageLockedAllocator<T>>;
    private:
        T cutoff;
        T squared_cutoff;
        T pot_shift;
        DeviceMDCell cell;
        DeviceCellList cellList;
        ForceBufferType forceBuffer;
        VirialBufferType virialBuffer;
        PageLockedVector swapBuffer;
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] __host__ __device__ T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] __host__ __device__ T force_functor(size_t i, size_t j, T r, T r2) const;

        [[nodiscard]] T potentialV(const MDCellType& hostCell) const;

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos);
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& hostCell);

        template<ExecutePolicy P>
        void forceAsync(const LatticeMatrix& lattice, const InvLatticeMatrix& invLattice, const PositionMatrix& cartesianPos, Vector auto& result);
        template<Vector V, ExecutePolicy P>
        void forceAsyncImpl(const LatticeMatrix& lattice, const InvLatticeMatrix& invLattice, const PositionMatrix& cartesianPos, V& result);
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result);
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell);
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos);
        [[nodiscard]] LatticeMatrix virial(const MDCellType& hostCell);
        void swap(This& __restrict obj) noexcept;

        __device__ void forceKernel();
        __device__ void postForceKernel();
        __device__ void virialKernel();
        __device__ void postVirialKernel();
        /* Getters */
        [[nodiscard]] __host__ __device__ const T& getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] __host__ __device__ const T& getSquaredCutoff() const noexcept { return squared_cutoff; }
        [[nodiscard]] __device__ const DeviceMDCell& getCell() const noexcept { return cell; }
        /* Setters */
        void setCutoff(T cutoff_);
    protected:
        device_obj() = default;
        device_obj(size_t numParticle);
        device_obj(size_t numParticle, T cutoff_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        size_t preParallel(
                const LatticeMatrix& lattice,
                const InvLatticeMatrix& invLattice,
                const PositionMatrix& cartesianPos,
                dim3& gridDims);
        [[nodiscard]] __device__ size_t forPairInCutoff(std::invocable<int, int, DeviceVector3D, T, T> auto fn) const;
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<PairModel<T>>> : public Traits<PairModel<T>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "PairModelImpl.cuh"
