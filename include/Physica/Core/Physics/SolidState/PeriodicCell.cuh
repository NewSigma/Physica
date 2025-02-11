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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh" // IWYU pragma: export
#include "PeriodicCell.h"

namespace Physica {
    template<Scalar T, unsigned int Dim>
    class device_obj<PeriodicCell<T, Dim>> : public Internal::PeriodicCellImpl {
        using host_obj = PeriodicCell<T, Dim>;
        using This = device_obj<host_obj>;
    public:
        using LatticeMatrix = device_obj<typename host_obj::LatticeMatrix>;
        using InvLatticeMatrix = device_obj<typename host_obj::InvLatticeMatrix>;
        using PositionMatrix = device_obj<typename host_obj::PositionMatrix>;
        using MomentumMatrix = PositionMatrix;
        using SearchRangeType = Array<ssize_t, Dim>;
    protected:
        using VectorType = device_obj<DenseVector<T, Dim>>;

        LatticeMatrix lattice;
        PositionMatrix pos;
        Type type;
    public:
        device_obj();
        device_obj(size_t numParticle, Type type_);
        device_obj(const host_obj& obj);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static unsigned int getDim() { return Dim; }
        [[nodiscard]] __host__ __device__ const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] __host__ __device__ const PositionMatrix& getPos() const noexcept { return pos; }
        [[nodiscard]] __host__ __device__ Type getType() const noexcept { return type; }
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return pos.getRow(); }
        [[nodiscard]] __device__ T getVolume() const noexcept { return getVolume(lattice); }
        /* Setters */
        void setLattice(LatticeMatrix new_lattice) { lattice = std::move(new_lattice); }
        /* Static members */
        [[nodiscard]] __device__ static T getVolume(const LatticeMatrix& lattice);
    };

    template<Scalar T, unsigned int Dim>
    device_obj<PeriodicCell<T, Dim>>::device_obj()
            : lattice(LatticeMatrix::unitMatrix(Dim))
            , pos()
            , type(Type::Direct) {}

    template<Scalar T, unsigned int Dim>
    device_obj<PeriodicCell<T, Dim>>::device_obj(size_t numParticle, Type type_)
            : device_obj() {
        pos.resize(numParticle, Dim);
        type = type_;
    }

    template<Scalar T, unsigned int Dim>
    device_obj<PeriodicCell<T, Dim>>::device_obj(const host_obj& obj)
            : lattice(obj.lattice)
            , pos(obj.pos)
            , type(obj.type) {}

    template<Scalar T, unsigned int Dim>
    void device_obj<PeriodicCell<T, Dim>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lattice.swap(obj.lattice);
        pos.swap(obj.pos);
        std::swap(type, obj.type);
    }

    template<Scalar T, unsigned int Dim>
    __device__ T device_obj<PeriodicCell<T, Dim>>::getVolume(const LatticeMatrix& lattice) {
        if constexpr (Dim == 1)
            return abs(lattice(0, 0));
        else if constexpr (Dim == 2)
            return (lattice.row(0).crossProduct(lattice.row(1))).compute().norm();
        else
            return abs(VectorType(lattice.row(0).crossProduct(lattice.row(1))) * lattice.row(2));
    }
}
