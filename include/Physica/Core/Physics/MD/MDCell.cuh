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

#include "MDCell.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.cuh"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim>
    class device_obj<MDCell<ScalarType, Dim>> : public device_obj<PeriodicCell<ScalarType, Dim>> {
        using host_obj = MDCell<ScalarType, Dim>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<PeriodicCell<ScalarType, Dim>>;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Type;
        using MassVector = device_obj<Vector<ScalarType>>;
    private:
        MassVector massVec;
        InvLatticeMatrix invLattice;
    public:
        device_obj() = default;
        explicit device_obj(size_t numParticle);
        device_obj(const host_obj& obj);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getDOF() const noexcept { return Dim * Base::getNumParticle(); }
        [[nodiscard]] __host__ __device__ const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] __host__ __device__ ScalarType getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] __host__ __device__ const InvLatticeMatrix& getInvLattice() const noexcept { return invLattice; }
        [[nodiscard]] __host__ __device__ constexpr static Type getType() noexcept { return Type::Cartesian; }
        /* Setters */
        void setLattice(const host_obj& cell);
    private:
        using Base::getType;
        using Base::setLattice;
    };

    template<class ScalarType, unsigned int Dim>
    device_obj<MDCell<ScalarType, Dim>>::device_obj(size_t numParticle)
            : Base(numParticle, Base::Type::Cartesian)
            , massVec(numParticle) {}

    template<class ScalarType, unsigned int Dim>
    device_obj<MDCell<ScalarType, Dim>>::device_obj(const host_obj& obj)
            : Base(obj)
            , massVec(obj.getMassVec())
            , invLattice(obj.getInvLattice()) {}

    template<class ScalarType, unsigned int Dim>
    void device_obj<MDCell<ScalarType, Dim>>::swap(device_obj& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        massVec.swap(obj.massVec);
        invLattice.swap(obj.invLattice);
    }

    template<class ScalarType, unsigned int Dim>
    void device_obj<MDCell<ScalarType, Dim>>::setLattice(const host_obj& cell) {
        Base::setLattice(cell.getLattice());
        invLattice = cell.getInvLattice();
    }
}
