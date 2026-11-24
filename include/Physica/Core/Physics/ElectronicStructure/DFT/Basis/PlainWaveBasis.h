/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/IndexVar/Index.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"
#include "Physica/Core/Scalar/Complex.h"

namespace Physica {
    template<Scalar T>
    class PlainWaveBasis {
        using ComplexType = Complex<T>;
        using GridType = DenseTensor<ComplexType, 3>;
        using LatticeMatrix = PeriodicCell<T, 3>::LatticeMatrix;

        GridType coeffGrid;
        LatticeMatrix repLatt;
        T cutEnergyPsi;
    public:
        PlainWaveBasis() = default;
        PlainWaveBasis(T cutEnergyPsi_, LatticeMatrix repLatt_);
        PlainWaveBasis(const PlainWaveBasis&) = default;
        PlainWaveBasis(PlainWaveBasis&&) noexcept = default;
        ~PlainWaveBasis() = default;
        /* Operators */
        PlainWaveBasis& operator=(PlainWaveBasis obj) noexcept;
        PlainWaveBasis& operator=(const Vector auto& statePsi);
        /* Operations */
        [[nodiscard]] Vector3D<T> makeWaveVector(size_t x, size_t y, size_t z) const noexcept;
        [[nodiscard]] Vector3D<T> makeWaveVector(Index3D index) const noexcept;
        void normalize();
        T calcNumElectron() const;

        template<RNG R>
        void random_normal();
        void swap(PlainWaveBasis& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getCoeffGrid(this auto&&) noexcept;
        [[nodiscard]] size_t getNumPlainWave() const noexcept { return coeffGrid.getSize(); }
        [[nodiscard]] Index3D getShape() const noexcept { return coeffGrid.getShape(); }
        [[nodiscard]] const LatticeMatrix& getRepLattice() const noexcept { return repLatt; }
        [[nodiscard]] T getCutEnergyPsi() const noexcept { return cutEnergyPsi; }
        /* Static members */
        static Index3D makeGridDim(T cutEnergy, const LatticeMatrix& repLatt);
        static Vector3D<T> makeWaveVector(const LatticeMatrix& repLatt, Index3D index, Index3D dim) noexcept;
        static void forKInBasis(const LatticeMatrix& repLatt, Index3D dim, std::invocable<Vector3D<T>> auto fn);
        static void forKIndexInBasis(const LatticeMatrix& repLatt, Index3D dim, std::invocable<Vector3D<T>, Index3D> auto fn);
    };

    template<Scalar T>
    PlainWaveBasis<T>::PlainWaveBasis(T cutEnergyPsi_, LatticeMatrix repLatt_)
            : repLatt(std::move(repLatt_))
            , cutEnergyPsi(cutEnergyPsi_) {
        coeffGrid.resize(makeGridDim(cutEnergyPsi, repLatt));
    }

    template<Scalar T>
    PlainWaveBasis<T>& PlainWaveBasis<T>::operator=(PlainWaveBasis<T> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    PlainWaveBasis<T>& PlainWaveBasis<T>::operator=(const Vector auto& statePsi) {
        assert(getNumPlainWave() == statePsi.getLength());
        auto flatGrid = coeffGrid.flatten();
        flatGrid = statePsi;
        return *this;
    }

    template<Scalar T>
    Vector3D<T> PlainWaveBasis<T>::makeWaveVector(size_t x, size_t y, size_t z) const noexcept {
        return makeWaveVector({x, y, z});
    }

    template<Scalar T>
    Vector3D<T> PlainWaveBasis<T>::makeWaveVector(Index3D index) const noexcept {
        return makeWaveVector(repLatt, index, getShape());
    }

    template<Scalar T>
    void PlainWaveBasis<T>::normalize() {
        auto vec = coeffGrid.flatten();
        vec *= sqrt(T(getNumPlainWave())) / vec.norm();
    }

    template<Scalar T>
    T PlainWaveBasis<T>::calcNumElectron() const {
        return coeffGrid.flatten().squaredNorm() / T(getNumPlainWave());
    }

    template<Scalar T>
    template<RNG R>
    void PlainWaveBasis<T>::random_normal() {
        coeffGrid.template random_normal<R>();
    }

    template<Scalar T>
    void PlainWaveBasis<T>::swap(PlainWaveBasis& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        coeffGrid.swap(obj.coeffGrid);
        repLatt.swap(obj.repLatt);
        cutEnergyPsi.swap(obj.cutEnergyPsi);
    }

    template<Scalar T>
    auto&& PlainWaveBasis<T>::getCoeffGrid(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).coeffGrid;
    }

    template<Scalar T>
    Index3D PlainWaveBasis<T>::makeGridDim(T cutEnergy, const LatticeMatrix& repLatt) {
        constexpr double factor = 2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck;
        const T maxWaveVec = sqrt(T(factor) * cutEnergy);
        const auto range = PeriodicCell<T, 3>::estimateRange(repLatt, maxWaveVec);
        return {size_t(2 * range[0] + 1), size_t(2 * range[1] + 1), size_t(2 * range[2] + 1)};
    }

    template<Scalar T>
    Vector3D<T> PlainWaveBasis<T>::makeWaveVector(const LatticeMatrix& repLatt, Index3D index, Index3D dim) noexcept {
        Vector3D<T> v{};
        for (int i = 0; i < 3; ++i) {
            const ssize_t index_i = index[i];
            const ssize_t dim_i = dim[i];
            assert(index_i < dim_i);
            v[i] = T(index_i > dim_i / 2 ? index_i - dim_i : index_i);
        }
        return repLatt.transpose() * v;
    }

    template<Scalar T>
    void PlainWaveBasis<T>::forKInBasis(const LatticeMatrix& repLatt, Index3D dim, std::invocable<Vector3D<T>> auto fn) {
        forND(dim, [&repLatt, dim, fn](Index3D index) {
            fn(makeWaveVector(repLatt, index, dim));
        });
    }

    template<Scalar T>
    void PlainWaveBasis<T>::forKIndexInBasis(const LatticeMatrix& repLatt, Index3D dim, std::invocable<Vector3D<T>, Index3D> auto fn) {
        forND(dim, [&repLatt, dim, fn](Index3D index) {
            fn(makeWaveVector(repLatt, index, dim), index);
        });
    }
}
