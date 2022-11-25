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
    class KSpaceGrid : private RSpaceGrid<typename std::conditional<T::isComplex, T, ComplexScalar<T>>::type> {
        constexpr static bool isComplex = T::isComplex;
        using RealType = typename T::RealType;
        using ComplexType = ComplexScalar<RealType>;
        using Base = RSpaceGrid<ComplexType>;
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
    public:
        using typename Base::Container;
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
        [[nodiscard]] ComplexType& operator()(ssize_t x, ssize_t y, ssize_t z);
        [[nodiscard]] const ComplexType& operator()(ssize_t x, ssize_t y, ssize_t z) const;
        [[nodiscard]] ComplexType& operator()(Index3D index) { return this->operator()(index[0], index[1], index[2]); }
        [[nodiscard]] const ComplexType& operator()(Index3D index) const { return this->operator()(index[0], index[1], index[2]); }
        /* Operations */
        [[nodiscard]] ComplexType calc(ssize_t x, ssize_t y, ssize_t z) const;
        void swap(KSpaceGrid& grid) noexcept { Base::swap(grid); }
        /* Getters */
        using Base::asVector;
        [[nodiscard]] inline size_t getSize() const noexcept;
        [[nodiscard]] ssize_t getDimX() const noexcept { return Base::getDimX() / 2U; }
        [[nodiscard]] ssize_t getDimY() const noexcept { return Base::getDimY() / 2U; }
        [[nodiscard]] ssize_t getDimZ() const noexcept { return isComplex ? (Base::getDimZ() / 2U) : (Base::getDimZ() - 1U); }
        [[nodiscard]] Index3D getDim() const noexcept { return {getDimX(), getDimY(), getDimZ()}; }
        /* Static members */
        static typename Base::Index3D makeGridDim(RealType cutEnergy, const LatticeMatrix& repCell);
        static KSpaceGrid<T> makeGrid(RealType cutEnergy, const LatticeMatrix& repCell);
        template<class Functor>
        inline static void forKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        static void forReducedKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        inline static void forKIndexInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
        template<class Functor>
        static void forReducedKIndexInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func);
        static void normalizeIndex(Index3D& index, const Index3D& range);
    protected:
        KSpaceGrid(Container data, size_t dimX_, size_t dimY_, size_t dimZ_);
    };

    template<class T>
    KSpaceGrid<T>::KSpaceGrid(size_t dimX, size_t dimY, size_t dimZ)
            : Base(2 * dimX + 1, 2 * dimY + 1, isComplex ? (2 * dimZ + 1) : (dimZ + 1)) {}

    template<class T>
    KSpaceGrid<T>::KSpaceGrid(Container data, size_t dimX, size_t dimY, size_t dimZ)
            : Base(std::move(data), 2 * dimX + 1, 2 * dimY + 1, isComplex ? (2 * dimZ + 1) : (dimZ + 1)) {}

    template<class T>
    KSpaceGrid<T>& KSpaceGrid<T>::operator=(KSpaceGrid<T> grid) noexcept {
        swap(grid);
        return *this;
    }

    template<class T>
    typename KSpaceGrid<T>::ComplexType& KSpaceGrid<T>::operator()(ssize_t x, ssize_t y, ssize_t z) {
        assert(-getDimX() <= x && x <= getDimX());
        assert(-getDimY() <= y && y <= getDimY());
        assert(-getDimZ() <= z && z <= getDimZ());
        if constexpr (isComplex)
            return Base::operator()(x + getDimX(), y + getDimY(), z + getDimZ());
        else {
            assert(z >= 0);
            return Base::operator()(x + getDimX(), y + getDimY(), z);
        }
    }

    template<class T>
    const typename KSpaceGrid<T>::ComplexType& KSpaceGrid<T>::operator()(ssize_t x, ssize_t y, ssize_t z) const {
        assert(-getDimX() <= x && x <= getDimX());
        assert(-getDimY() <= y && y <= getDimY());
        assert(-getDimZ() <= z && z <= getDimZ());
        if constexpr (isComplex)
            return Base::operator()(x + getDimX(), y + getDimY(), z + getDimZ());
        else {
            assert(z >= 0);
            return Base::operator()(x + getDimX(), y + getDimY(), z);
        }
    }

    template<class T>
    typename KSpaceGrid<T>::ComplexType KSpaceGrid<T>::calc(ssize_t x, ssize_t y, ssize_t z) const {
        const ComplexType temp = isComplex ? this->operator()(x, y, z) : this->operator()(x, y, std::abs(z));
        return z >= 0 ? temp : temp.conjugate();
    }

    template<class T>
    size_t KSpaceGrid<T>::getSize() const noexcept {
        if constexpr (isComplex)
            return Base::getDimX() * Base::getDimY() * Base::getDimZ();
        return Base::getDimX() * Base::getDimY() * (Base::getDimZ() * 2 - 1);
    }

    template<class T>
    typename KSpaceGrid<T>::Base::Index3D KSpaceGrid<T>::makeGridDim(RealType cutEnergy, const LatticeMatrix& repCell) {
        constexpr double factor = 2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck;
        const RealType maxWaveVec = sqrt(RealType(factor) * cutEnergy);
        const auto range = PeriodicCell<RealType, 3>::estimateRange(repCell, maxWaveVec);
        return {(size_t)range[0], (size_t)range[1], (size_t)range[2]};
    }

    template<class T>
    KSpaceGrid<T> KSpaceGrid<T>::makeGrid(RealType cutEnergy, const LatticeMatrix& repCell) {
        const auto range = makeGridDim(cutEnergy, repCell);
        return KSpaceGrid<T>(range[0], range[1], range[2]);
    }

    template<class T>
    template<class Functor>
    inline void KSpaceGrid<T>::forKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func) {
        PeriodicCell<RealType, 3>::forCellInRange(dim, repLatt, func);
    }

    template<class T>
    template<class Functor>
    void KSpaceGrid<T>::forReducedKInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func) {
        if constexpr (isComplex)
            forKInGrid(dim, repLatt, func);
        else {
            auto a1 = repLatt.row(0);
            auto a2 = repLatt.row(1);
            auto a3 = repLatt.row(2);

            VectorType v1, v2, v3;
            for (ssize_t x = -dim[0]; x <= dim[0]; ++x) {
                v1 = RealType(x) * a1.asVector();
                for (ssize_t y = -dim[1]; y <= dim[1]; ++y) {
                    v2 = v1 + RealType(y) * a2.asVector();
                    for (ssize_t z = 0; z <= dim[2]; ++z) {
                        v3 = v2 + RealType(z) * a3.asVector();
                        func(v3);
                    }
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
            v1 = RealType(x) * a1.asVector();
            for (ssize_t y = -dim[1]; y <= dim[1]; ++y) {
                v2 = v1 + RealType(y) * a2.asVector();
                for (ssize_t z = -dim[2]; z <= dim[2]; ++z) {
                    v3 = v2 + RealType(z) * a3.asVector();
                    func(v3, Index3D{x, y, z});
                }
            }
        }
    }

    template<class T>
    template<class Functor>
    void KSpaceGrid<T>::forReducedKIndexInGrid(Index3D dim, const LatticeMatrix& repLatt, Functor func) {
        if constexpr (isComplex)
            forKIndexInGrid(dim, repLatt, func);
        else {
            auto a1 = repLatt.row(0);
            auto a2 = repLatt.row(1);
            auto a3 = repLatt.row(2);

            VectorType v1, v2, v3;
            for (ssize_t x = -dim[0]; x <= dim[0]; ++x) {
                v1 = RealType(x) * a1.asVector();
                for (ssize_t y = -dim[1]; y <= dim[1]; ++y) {
                    v2 = v1 + RealType(y) * a2.asVector();
                    for (ssize_t z = 0; z <= dim[2]; ++z) {
                        v3 = v2 + RealType(z) * a3.asVector();
                        func(v3, Index3D{x, y, z});
                    }
                }
            }
        }
    }

    template<class T>
    void KSpaceGrid<T>::normalizeIndex(Index3D& index, const Index3D& range) {
        for (int i = 0; i < 3; ++i) {
            if (index[i] > range[i])
                index[i] -= range[i];
            else if (index[i] < -range[i])
                index[i] += range[i];
            assert(0 <= index[i] && index[1] < range[i]);
        }
    }
}
