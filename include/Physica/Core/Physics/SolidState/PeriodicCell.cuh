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

#include "PeriodicCell.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim>
    class device_obj<PeriodicCell<ScalarType, Dim>> : public Internal::PeriodicCellImpl {
        using host_obj = PeriodicCell<ScalarType, Dim>;
        using This = device_obj<host_obj>;
    public:
        using LatticeMatrix = device_obj<typename host_obj::LatticeMatrix>;
        using PositionMatrix = device_obj<typename host_obj::PositionMatrix>;
        using MomentumMatrix = PositionMatrix;
        using SearchRangeType = Utils::Array<ssize_t, Dim>;
    protected:
        using InvLatticeMatrix = device_obj<typename host_obj::InvLatticeMatrix>;
        using VectorType = device_obj<Vector<ScalarType, Dim>>;

        LatticeMatrix lattice;
        PositionMatrix pos;
        Type type;
    public:
        device_obj();
        device_obj(size_t numParticle, Type type_);
        device_obj(const host_obj& obj);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static unsigned int getDim() { return Dim; }
        [[nodiscard]] __host__ __device__ const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] __host__ __device__ const PositionMatrix& getPos() const noexcept { return pos; }
        [[nodiscard]] __host__ __device__ Type getType() const noexcept { return type; }
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return pos.getRow(); }
        /* Setters */
        void setLattice(LatticeMatrix new_lattice) { lattice = new_lattice; }
    };

    template<class ScalarType, unsigned int Dim>
    device_obj<PeriodicCell<ScalarType, Dim>>::device_obj()
            : lattice(LatticeMatrix::unitMatrix(Dim))
            , pos()
            , type(Type::Direct) {}

    template<class ScalarType, unsigned int Dim>
    device_obj<PeriodicCell<ScalarType, Dim>>::device_obj(size_t numParticle, Type type_)
            : device_obj() {
        pos.resize(numParticle, Dim);
        type = type_;
    }

    template<class ScalarType, unsigned int Dim>
    device_obj<PeriodicCell<ScalarType, Dim>>::device_obj(const host_obj& obj)
            : lattice(obj.lattice)
            , pos(obj.pos)
            , type(obj.type) {}

    template<class ScalarType, unsigned int Dim>
    void device_obj<PeriodicCell<ScalarType, Dim>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lattice.swap(obj.lattice);
        pos.swap(obj.pos);
        std::swap(type, obj.type);
    }
}
