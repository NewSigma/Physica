/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/KSpaceGrid.h"

namespace Physica::Core {
    /**
     * \class PWBaseWave
     * Plain wave base wave function
     */
    template<class ScalarType>
    class PWBaseWave : public KSpaceGrid<RSpaceGrid<ComplexScalar<ScalarType>>> {
    public:
        using ComplexType = ComplexScalar<ScalarType>;
        using Base = KSpaceGrid<RSpaceGrid<ComplexType>>;
    private:
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;

        LatticeMatrix repCell;
    public:
        PWBaseWave() = default;
        PWBaseWave(ScalarType cutEnergy, LatticeMatrix repCell_);
        PWBaseWave(const PWBaseWave&) = default;
        PWBaseWave(PWBaseWave&&) noexcept = default;
        ~PWBaseWave() = default;
        /* Operators */
        using Base::operator();
        PWBaseWave& operator=(PWBaseWave wave) noexcept;
        [[nodiscard]] ComplexType operator()(Vector<ScalarType, 3> k, Vector<ScalarType, 3> r) const;
        template<class VectorType>
        PWBaseWave& operator=(const RValueVector<VectorType>& newCoeffs);
        template<class T>
        friend std::ostream& operator<<(std::ostream& os, const PWBaseWave<T>& wave);
        template<class T>
        friend std::istream& operator>>(std::istream& is, PWBaseWave<T>& wave);
        /* Operations */
        void swap(PWBaseWave& wave) noexcept;
        /* Getters */
        using Base::flatten;
        using Base::getDimX;
        using Base::getDimY;
        using Base::getDimZ;
        using Base::getDim;
        [[nodiscard]] size_t getNumPlaneWave() const noexcept { return Base::getSize(); }
        [[nodiscard]] Vector<ScalarType, 3> getWaveVector(size_t x, size_t y, size_t z) const noexcept;
        [[nodiscard]] ScalarType getKinetic(size_t x, size_t y, size_t z) const noexcept;
    };

    template<class ScalarType>
    PWBaseWave<ScalarType>::PWBaseWave(ScalarType cutEnergy, LatticeMatrix repCell_) : repCell(std::move(repCell_)) {
        Base::operator=(Base(KSpaceGrid<ComplexType>::makeGridDim(cutEnergy, repCell)));
    }

    template<class ScalarType>
    PWBaseWave<ScalarType>& PWBaseWave<ScalarType>::operator=(PWBaseWave<ScalarType> wave) noexcept {
        swap(wave);
        return *this;
    }

    template<class ScalarType>
    typename PWBaseWave<ScalarType>::ComplexType PWBaseWave<ScalarType>::operator()(
            Vector<ScalarType, 3> k, Vector<ScalarType, 3> r) const {
        using Index3D = typename Base::Index3D;
        ComplexType result = ComplexType(0);
        Base::forKIndexInGrid(*this, repCell, [this, &result, &k, &r](Vector<ScalarType, 3> K, Index3D index) {
            const ScalarType phase = (k + K) * r;
            ScalarType s, c;
            sincos(phase, s, c);
            result += this->operator()(index) * ComplexType(c, s);
        });
        return result;
    }

    template<class ScalarType>
    template<class VectorType>
    PWBaseWave<ScalarType>& PWBaseWave<ScalarType>::operator=(const RValueVector<VectorType>& newCoeffs) {
        flatten() = newCoeffs;
        return *this;
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const PWBaseWave<T>& wave) {
        os << wave.repCell;
        os << static_cast<const typename PWBaseWave<T>::Base&>(wave);
        return os;
    }

    template<class T>
    std::istream& operator>>(std::istream& is, PWBaseWave<T>& wave) {
        is >> wave.repCell;
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        is >> static_cast<typename PWBaseWave<T>::Base&>(wave);
        return is;
    }

    template<class ScalarType>
    void PWBaseWave<ScalarType>::swap(PWBaseWave& wave) noexcept {
        Base::swap(wave);
        repCell.swap(wave.repCell);
    }

    template<class ScalarType>
    Vector<ScalarType, 3> PWBaseWave<ScalarType>::getWaveVector(size_t x, size_t y, size_t z) const noexcept {
        const ssize_t x1 = x > getDimX() / 2 ? -ssize_t(getDimX() - x) : x;
        const ssize_t y1 = y > getDimY() / 2 ? -ssize_t(getDimY() - y) : y;
        const ssize_t z1 = z > getDimZ() / 2 ? -ssize_t(getDimZ() - z) : z;
        return repCell.row(0).asVector() * ScalarType(x1) +
               repCell.row(1).asVector() * ScalarType(y1) +
               repCell.row(2).asVector() * ScalarType(z1);
    }

    template<class ScalarType>
    ScalarType PWBaseWave<ScalarType>::getKinetic(size_t x, size_t y, size_t z) const noexcept {
        assert(x < getDimX());
        assert(y < getDimY());
        assert(z < getDimZ());
        constexpr double factor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck / PhyConst<AU>::electronMass * 0.5;
        return getWaveVector(x, y, z).squaredNorm() * ScalarType(factor);
    }
}
