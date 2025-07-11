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

#include "EmptyForceModel.h"

namespace Physica {
    template<Scalar T, unsigned int Dim>
    class Hamonic : private EmptyForceModel<T, Dim> {
        using Base = EmptyForceModel<T, Dim>;
        using This = Hamonic<T, Dim>;
    public:
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::ForceConstMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
    private:
        PositionMatrix sites;
        VectorND<T> springCoeffs;
    public:
        Hamonic(PositionMatrix sites_, VectorND<T> springCoeffs_);
        Hamonic(const This&) = default;
        Hamonic(This&&) noexcept = default;
        ~Hamonic() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T potentialV(const MDCellType& cell) const;

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result) const;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<P>(cell); }
        using Base::force_long;

        [[nodiscard]] T forceConst([[maybe_unused]] const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return sites.getRow(); }
        [[nodiscard]] const VectorND<T>& getSpringCoeffs() const noexcept { return springCoeffs; }
    };

    template<Scalar T, unsigned int Dim>
    Hamonic<T, Dim>::Hamonic(PositionMatrix sites_, VectorND<T> springCoeffs_)
            : sites(std::move(sites_))
            , springCoeffs(std::move(springCoeffs_)) {
        assert(springCoeffs.getLength() == getNumParticle() && "[Error]: Number of particles is not consistent");
    }

    template<Scalar T, unsigned int Dim>
    T Hamonic<T, Dim>::potentialV(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        T result = 0;
        for (size_t i = 0; i < getNumParticle(); ++i)
            result += springCoeffs[i] * cell.minDistVector(sites.row(i), i).squaredNorm();
        return result * T(0.5);
    }

    template<Scalar T, unsigned int Dim>
    template<ExecutePolicy P>
    VectorND<T> Hamonic<T, Dim>::force(const MDCellType& cell) const {
        VectorND<T> result(cell.getDOF());
        forceAsync<P>(cell, result);
        return result;
    }

    template<Scalar T, unsigned int Dim>
    template<ExecutePolicy P>
    void Hamonic<T, Dim>::forceAsync(const MDCellType& cell, Vector auto& result) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        result = T(0);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            auto force_i = result.template segment<Dim>(i * Dim, (i + 1) * Dim);
            force_i = -springCoeffs[i] * cell.minDistVector(sites.row(i), i);
        }
    }

    template<Scalar T, unsigned int Dim>
    T Hamonic<T, Dim>::forceConst([[maybe_unused]] const MDCellType& cell, size_t dof1, size_t dof2) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        assert(dof1 < cell.getDOF() && dof2 < cell.getDOF() && "[Error]: Index overflow");
        if (dof1 == dof2) {
            const size_t atom = dof1 / Dim;
            return springCoeffs[atom];
        }
        return T(0);
    }

    template<Scalar T, unsigned int Dim>
    auto Hamonic<T, Dim>::forceConst(const MDCellType& cell) const -> ForceConstMatrix {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        const size_t dof = cell.getDOF();
        ForceConstMatrix result(dof, T(0));
        for (size_t i = 0; i < dof; ++i)
            result(i, i) = forceConst(cell, i, i);
        return result;
    }

    template<Scalar T, unsigned int Dim>
    auto Hamonic<T, Dim>::virial(const MDCellType& cell) const -> LatticeMatrix {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const DenseVector<T, Dim> delta = cell.minDistVector(sites.row(i), i);
            const DenseVector<T, Dim> temp = springCoeffs[i] * delta;
            result += temp * delta.transpose();
        }
        result *= -reciprocal(cell.getVolume());
        return result;
    }

    template<Scalar T, unsigned int Dim>
    void Hamonic<T, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        sites.swap(obj.sites);
        springCoeffs.swap(obj.springCoeffs);
    }
}

namespace Physica {
    template<Scalar T, unsigned int Dim>
    class Traits<Hamonic<T, Dim>> {
    public:
        constexpr static bool IsPeriodBoundary = true;
        constexpr static bool IsContractable = false;
    };
}
