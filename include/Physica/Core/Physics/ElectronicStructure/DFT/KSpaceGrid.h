/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Physics/PeriodicCell.h"
#include "RSpaceGrid.h"

namespace Physica::Core {
    template<class T>
    class KSpaceGrid : private RSpaceGrid<T> {
        using Base = RSpaceGrid<T>;
        using ScalarType = typename T::RealType;
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using typename Base::Container;
    public:
        using Index3D = Utils::Array<ssize_t, 3>;
        using VectorType = Vector<typename LatticeMatrix::ScalarType, 3>;
    public:
        KSpaceGrid() = default;
        KSpaceGrid(size_t dimX, size_t dimY, size_t dimZ);
        KSpaceGrid(const KSpaceGrid&) = default;
        KSpaceGrid(KSpaceGrid&&) noexcept = default;
        ~KSpaceGrid() = default;
        /* Operators */
        KSpaceGrid& operator=(KSpaceGrid grid) noexcept;
        [[nodiscard]] T& operator()(ssize_t x, ssize_t y, ssize_t z);
        [[nodiscard]] const T& operator()(ssize_t x, ssize_t y, ssize_t z) const;
        [[nodiscard]] T& operator()(Index3D index) { return this->operator()(index[0], index[1], index[2]); }
        [[nodiscard]] const T& operator()(Index3D index) const { return this->operator()(index[0], index[1], index[2]); }
        /* Operations */
        void swap(KSpaceGrid& grid) noexcept { Base::swap(grid); }
        /* Getters */
        using Base::asVector;
        [[nodiscard]] size_t getSize() const noexcept { return Base::getDimX() * Base::getDimY() * (Base::getDimZ() * 2 - 1); }
        [[nodiscard]] ssize_t getDimX() const noexcept { return Base::getDimX() / 2U; }
        [[nodiscard]] ssize_t getDimY() const noexcept { return Base::getDimY() / 2U; }
        [[nodiscard]] ssize_t getDimZ() const noexcept { return Base::getDimZ() - 1U; }
        [[nodiscard]] Index3D getDim() const noexcept { return {getDimX(), getDimY(), getDimZ()}; }
        /* Static members */
        static typename Base::Index3D makeGridDim(ScalarType cutEnergy, const LatticeMatrix& repCell);
        static KSpaceGrid<T> makeGrid(ScalarType cutEnergy, const LatticeMatrix& repCell);
        template<class Functor>
        inline static void forKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        static void forReducedKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        inline static void forKIndexInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
    protected:
        KSpaceGrid(Container data, size_t dimX_, size_t dimY_, size_t dimZ_);
    };

    template<class T>
    KSpaceGrid<T>::KSpaceGrid(size_t dimX, size_t dimY, size_t dimZ) : Base(2 * dimX + 1, 2 * dimY + 1, dimZ + 1) {}

    template<class T>
    KSpaceGrid<T>::KSpaceGrid(Container data, size_t dimX, size_t dimY, size_t dimZ)
            : Base(std::move(data), 2 * dimX + 1, 2 * dimY + 1, dimZ + 1) {}

    template<class T>
    KSpaceGrid<T>& KSpaceGrid<T>::operator=(KSpaceGrid<T> grid) noexcept {
        swap(grid);
        return *this;
    }

    template<class T>
    T& KSpaceGrid<T>::operator()(ssize_t x, ssize_t y, ssize_t z) {
        assert(-getDimX() <= x && x <= getDimX());
        assert(-getDimY() <= y && y <= getDimY());
        assert(-getDimZ() <= z && z <= getDimZ());
        return Base::operator()(x + getDimX(), y + getDimY(), std::abs(z));
    }

    template<class T>
    const T& KSpaceGrid<T>::operator()(ssize_t x, ssize_t y, ssize_t z) const {
        assert(-getDimX() <= x && x <= getDimX());
        assert(-getDimY() <= y && y <= getDimY());
        assert(-getDimZ() <= z && z <= getDimZ());
        return Base::operator()(x + getDimX(), y + getDimY(), std::abs(z));
    }

    template<class T>
    typename RSpaceGrid<T>::Index3D KSpaceGrid<T>::makeGridDim(ScalarType cutEnergy, const LatticeMatrix& repCell) {
        constexpr double factor = 2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck;
        const ScalarType maxWaveVec = sqrt(ScalarType(factor) * cutEnergy);
        const auto range = PeriodicCell<ScalarType, 3>::estimateRange(repCell, maxWaveVec);
        return {(size_t)range[0], (size_t)range[1], (size_t)range[2]};
    }

    template<class T>
    KSpaceGrid<T> KSpaceGrid<T>::makeGrid(ScalarType cutEnergy, const LatticeMatrix& repCell) {
        const auto range = makeGridDim(cutEnergy, repCell);
        return KSpaceGrid<T>(range[0], range[1], range[2]);
    }

    template<class T>
    template<class Functor>
    inline void KSpaceGrid<T>::forKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func) {
        PeriodicCell<ScalarType, 3>::forCellInRange(dim, repLatt, func);
    }

    template<class T>
    template<class Functor>
    void KSpaceGrid<T>::forReducedKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func) {
        auto a1 = repLatt.row(0);
        auto a2 = repLatt.row(1);
        auto a3 = repLatt.row(2);

        VectorType v1, v2, v3;
        for (ssize_t x = -dim[0]; x <= dim[0]; ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (ssize_t y = -dim[1]; y <= dim[1]; ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (ssize_t z = 0; z <= dim[2]; ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3);
                }
            }
        }
    }

    template<class T>
    template<class Functor>
    void KSpaceGrid<T>::forKIndexInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func) {
        auto a1 = repLatt.row(0);
        auto a2 = repLatt.row(1);
        auto a3 = repLatt.row(2);

        VectorType v1, v2, v3;
        for (ssize_t x = -dim[0]; x <= dim[0]; ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (ssize_t y = -dim[1]; y <= dim[1]; ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (ssize_t z = -dim[2]; z <= dim[2]; ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3, Index3D{x, y, z});
                }
            }
        }
    }
}
