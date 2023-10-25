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

#include "PairModel.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Physics/MD/MDCell.cuh"
#include "Physica/Core/Physics/MD/MDImpl/CellList.cuh"
#include "Physica/Utils/Allocator/PageLockedAllocator.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<PairModel<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        using host_obj = PairModel<Derived>;
        using This = device_obj<PairModel<Derived>>;
        using Base = Utils::CRTPBase<device_obj<Derived>>;
        using TraitType = Internal::Traits<Derived>;

        constexpr static int Dim = 3;
        constexpr static bool IsPotDependOnAtomIndex = TraitType::IsPotDependOnAtomIndex;
    public:
        using ScalarType = typename TraitType::ScalarType;
        using MDCellType = MDCell<ScalarType>;
        using DeviceMDCell = device_obj<MDCellType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using CellListType = CellList<ScalarType>;
        using DeviceCellList = device_obj<CellListType>;
        using Index3D = typename GridBase::Index3D;
        using Vector3D = device_obj<Vector<ScalarType, 3>>;
        using DeviceMatrix = device_obj<DenseMatrix<ScalarType>>;
        using PageLockedVector = Vector<ScalarType, Dynamic, Dynamic, Utils::PageLockedAllocator<ScalarType>>;
    private:
        ScalarType cutoff;
        ScalarType squared_cutoff;
        ScalarType pot_shift;
        DeviceMDCell cell;
        DeviceCellList cellList;
        DeviceMatrix forceBuffer;
        PageLockedVector swapBuffer;
    public:
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] __host__ __device__ inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] __host__ __device__ inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;

        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& hostCell) const;

        template<class Executor, bool IsSmallCell> [[nodiscard]] Vector<ScalarType> force(const MDCellType& hostCell);
        template<class VectorType, class Executor, bool IsSmallCell>
        void forceAsync(const MDCellType& hostCell, ContinuousVector<VectorType>& result);
        template<class VectorType, class Executor>
        inline void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result);
        template<class Executor>
        [[nodiscard]] inline Vector<ScalarType> force_short(const MDCellType& cell);
        template<class Executor>
        [[nodiscard]] inline Vector<ScalarType> force_long(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& hostCell) const;
        void swap(device_obj& obj) noexcept;

        template<bool IsSmallCell>
        __device__ void forceKernelImpl();
        __device__ void postForceKernelImpl();
        /* Getters */
        [[nodiscard]] __device__ const DeviceMDCell& getCell() const noexcept { return cell; }
        [[nodiscard]] bool isSmallCell(const MDCellType& cell) const noexcept { return Base::getDerived().isSmallCell(cell); }
    protected:
        device_obj() = default;
        device_obj(size_t numParticle, ScalarType cutoff_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
    };
}

#include "PairModelImpl.cuh"
