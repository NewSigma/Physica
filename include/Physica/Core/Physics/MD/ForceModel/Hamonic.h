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

#include "EmptyForceModel.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim> class Hamonic;

    namespace Internal {
        template<class ScalarType, unsigned int Dim>
        class Traits<Hamonic<ScalarType, Dim>> {
        public:
            constexpr static bool IsPeriodBoundary = true;
        };
    }

    template<class ScalarType, unsigned int Dim>
    class Hamonic : private EmptyForceModel<ScalarType, Dim> {
        using Base = EmptyForceModel<ScalarType, Dim>;
    public:
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::ForceConstMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
    private:
        PositionMatrix sites;
        Vector<ScalarType> springCoeffs;
    public:
        Hamonic(PositionMatrix sites_, Vector<ScalarType> springCoeffs_);
        Hamonic(const Hamonic&) = default;
        Hamonic(Hamonic&&) noexcept = default;
        ~Hamonic() = default;
        /* Operators */
        Hamonic& operator=(Hamonic obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;

        template<class Executor, bool IsSmallCell = false>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor, bool IsSmallCell = false>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        using Base::force_long;

        [[nodiscard]] ScalarType forceConst([[maybe_unused]] const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(Hamonic& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return sites.getRow(); }
    };

    template<class ScalarType, unsigned int Dim>
    Hamonic<ScalarType, Dim>::Hamonic(PositionMatrix sites_, Vector<ScalarType> springCoeffs_)
            : sites(std::move(sites_))
            , springCoeffs(std::move(springCoeffs_)) {
        assert(springCoeffs.getLength() == getNumParticle() && "[Error]: Number of particles is not consistent");
    }

    template<class ScalarType, unsigned int Dim>
    ScalarType Hamonic<ScalarType, Dim>::potentialEnergy(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        ScalarType result = 0;
        for (size_t i = 0; i < getNumParticle(); ++i)
            result += springCoeffs[i] * cell.minDistVector(sites.row(i), i).squaredNorm();
        return result * ScalarType(0.5);
    }

    template<class ScalarType, unsigned int Dim>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> Hamonic<ScalarType, Dim>::force(const MDCellType& cell) const {
        Vector<ScalarType> result(cell.getDOF());
        forceAsync<Vector<ScalarType>, Executor, IsSmallCell>(cell, result);
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    template<class VectorType, class Executor, bool IsSmallCell>
    void Hamonic<ScalarType, Dim>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        result = ScalarType(0);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            auto force_i = result.template segment<Dim>(i * Dim, (i + 1) * Dim);
            force_i = -springCoeffs[i] * cell.minDistVector(sites.row(i), i);
        }
    }

    template<class ScalarType, unsigned int Dim>
    ScalarType Hamonic<ScalarType, Dim>::forceConst([[maybe_unused]] const MDCellType& cell, size_t dof1, size_t dof2) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        assert(dof1 < cell.getDOF() && dof2 < cell.getDOF() && "[Error]: Index overflow");
        if (dof1 == dof2) {
            const size_t atom = dof1 / Dim;
            return springCoeffs[atom];
        }
        return ScalarType(0);
    }

    template<class ScalarType, unsigned int Dim>
    typename Hamonic<ScalarType, Dim>::ForceConstMatrix
    Hamonic<ScalarType, Dim>::forceConst(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        const size_t dof = cell.getDOF();
        ForceConstMatrix result(dof, ScalarType(0));
        for (size_t i = 0; i < dof; ++i)
            result(i, i) = forceConst(cell, i, i);
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    typename Hamonic<ScalarType, Dim>::LatticeMatrix Hamonic<ScalarType, Dim>::virial(const MDCellType& cell) const {
        assert(cell.getNumParticle() == getNumParticle() && "[Error]: Number of particles is not consistent");
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Vector<ScalarType, Dim> delta = cell.minDistVector(sites.row(i), i);
            const Vector<ScalarType, Dim> temp = springCoeffs[i] * delta;
            result += temp * delta.transpose();
        }
        result *= -reciprocal(cell.getVolume());
        return result;
    }

    template<class ScalarType, unsigned int Dim>
    void Hamonic<ScalarType, Dim>::swap(Hamonic& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        sites.swap(obj.sites);
        springCoeffs.swap(obj.springCoeffs);
    }
}
