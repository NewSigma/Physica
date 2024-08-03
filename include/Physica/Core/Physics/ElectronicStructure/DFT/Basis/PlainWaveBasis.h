/*
 * Copyright 2021-2023 Weibo He.
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

#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"

namespace Physica::Core {
    template<class ScalarType>
    class PlainWaveBasis {
        using ComplexType = ComplexScalar<ScalarType>;
        using GridType = RSpaceGrid<ComplexType>;
        using LatticeMatrix = typename PeriodicCell<ScalarType, 3>::LatticeMatrix;
        using Index3D = typename GridBase::Index3D;
        using Vector3D = Vector<ScalarType, 3>;

        GridType coeffGrid;
        LatticeMatrix repLatt;
        ScalarType cutEnergyPsi;
    public:
        PlainWaveBasis() = default;
        PlainWaveBasis(ScalarType cutEnergyPsi_, LatticeMatrix repLatt_);
        PlainWaveBasis(const PlainWaveBasis&) = default;
        PlainWaveBasis(PlainWaveBasis&&) noexcept = default;
        ~PlainWaveBasis() = default;
        /* Operators */
        PlainWaveBasis& operator=(PlainWaveBasis obj) noexcept;
        template<class VectorType>
        PlainWaveBasis& operator=(const RValueVector<VectorType>& statePsi);
        /* Operations */
        [[nodiscard]] inline Vector3D makeWaveVector(size_t x, size_t y, size_t z) const noexcept;
        [[nodiscard]] inline Vector3D makeWaveVector(Index3D index) const noexcept;
        inline void normalize();
        ScalarType calcNumElectron() const;

        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        void swap(PlainWaveBasis& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const GridType& getCoeffGrid() const noexcept { return coeffGrid; }
        [[nodiscard]] size_t getNumPlainWave() const noexcept { return coeffGrid.getSize(); }
        [[nodiscard]] Index3D getDim() const noexcept { return coeffGrid.getDim(); }
        [[nodiscard]] const LatticeMatrix& getRepLattice() const noexcept { return repLatt; }
        [[nodiscard]] ScalarType getCutEnergyPsi() const noexcept { return cutEnergyPsi; }
        /* Static members */
        static Index3D makeGridDim(ScalarType cutEnergy, const LatticeMatrix& repLatt);
        static Vector3D makeWaveVector(const LatticeMatrix& repLatt, Index3D index, Index3D dim) noexcept;
        template<class Functor>
        inline static void forKInBasis(const LatticeMatrix& repLatt, Index3D dim, Functor func);
        template<class Functor>
        inline static void forKIndexInBasis(const LatticeMatrix& repLatt, Index3D dim, Functor func);
    };

    template<class ScalarType>
    PlainWaveBasis<ScalarType>::PlainWaveBasis(ScalarType cutEnergyPsi_, LatticeMatrix repLatt_)
            : repLatt(std::move(repLatt_))
            , cutEnergyPsi(cutEnergyPsi_) {
        coeffGrid.resize(makeGridDim(cutEnergyPsi, repLatt));
    }

    template<class ScalarType>
    PlainWaveBasis<ScalarType>& PlainWaveBasis<ScalarType>::operator=(PlainWaveBasis<ScalarType> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    template<class VectorType>
    PlainWaveBasis<ScalarType>& PlainWaveBasis<ScalarType>::operator=(const RValueVector<VectorType>& statePsi) {
        assert(getNumPlainWave() == statePsi.getLength());
        auto flatGrid = coeffGrid.flatten();
        flatGrid = statePsi;
        return *this;
    }

    template<class ScalarType>
    inline typename PlainWaveBasis<ScalarType>::Vector3D
    PlainWaveBasis<ScalarType>::makeWaveVector(size_t x, size_t y, size_t z) const noexcept {
        return makeWaveVector({x, y, z});
    }

    template<class ScalarType>
    inline typename PlainWaveBasis<ScalarType>::Vector3D
    PlainWaveBasis<ScalarType>::makeWaveVector(Index3D index) const noexcept {
        return makeWaveVector(repLatt, index, getDim());
    }

    template<class ScalarType>
    inline void PlainWaveBasis<ScalarType>::normalize() {
        auto vec = coeffGrid.flatten();
        vec *= sqrt(ScalarType(getNumPlainWave())) / vec.norm();
    }

    template<class ScalarType>
    ScalarType PlainWaveBasis<ScalarType>::calcNumElectron() const {
        return coeffGrid.flatten().squaredNorm() / ScalarType(getNumPlainWave());
    }

    template<class ScalarType>
    template<class RandomGenerator>
    inline void PlainWaveBasis<ScalarType>::random_normal(RandomGenerator& gen) {
        coeffGrid.template random_normal<RandomGenerator>(gen);
    }

    template<class ScalarType>
    void PlainWaveBasis<ScalarType>::swap(PlainWaveBasis& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        coeffGrid.swap(obj.coeffGrid);
        repLatt.swap(obj.repLatt);
        cutEnergyPsi.swap(obj.cutEnergyPsi);
    }

    template<class ScalarType>
    typename GridBase::Index3D PlainWaveBasis<ScalarType>::makeGridDim(ScalarType cutEnergy, const LatticeMatrix& repLatt) {
        constexpr double factor = 2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck;
        const ScalarType maxWaveVec = sqrt(ScalarType(factor) * cutEnergy);
        const auto range = PeriodicCell<ScalarType, 3>::estimateRange(repLatt, maxWaveVec);
        return {size_t(2 * range[0] + 1), size_t(2 * range[1] + 1), size_t(2 * range[2] + 1)};
    }

    template<class ScalarType>
    typename PlainWaveBasis<ScalarType>::Vector3D
    PlainWaveBasis<ScalarType>::makeWaveVector(const LatticeMatrix& repLatt, Index3D index, Index3D dim) noexcept {
        Vector3D v{};
        for (int i = 0; i < 3; ++i) {
            const ssize_t index_i = index[i];
            const ssize_t dim_i = dim[i];
            assert(index_i < dim_i);
            v[i] = ScalarType(index_i > dim_i / 2 ? index_i - dim_i : index_i);
        }
        return repLatt.transpose() * v;
    }

    template<class ScalarType>
    template<class Functor>
    inline void PlainWaveBasis<ScalarType>::forKInBasis(const LatticeMatrix& repLatt, Index3D dim, Functor func) {
        GridBase::forIndexInGrid(dim, [&repLatt, dim, func](Index3D index) {
            func(makeWaveVector(repLatt, index, dim));
        });
    }

    template<class ScalarType>
    template<class Functor>
    inline void PlainWaveBasis<ScalarType>::forKIndexInBasis(const LatticeMatrix& repLatt, Index3D dim, Functor func) {
        GridBase::forIndexInGrid(dim, [&repLatt, dim, func](Index3D index) {
            func(makeWaveVector(repLatt, index, dim), index);
        });
    }
}
