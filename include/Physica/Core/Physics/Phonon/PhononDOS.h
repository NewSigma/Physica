/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/Math/Calculus/PDE/FEM/Element/CuboidLinear.h"
#include "Physica/Core/Math/Geometry/CubeCross.h"
#include "PhononSolver.h"

namespace Physica::Core {
    template<class ScalarType>
    class PhononDOS {
    public:
        using SolverType = PhononSolver<ScalarType>;
        using ElementType = CuboidLinear<ScalarType>;
        using Index3D = typename SolverType::Index3D;
        using Vector3D = typename SolverType::Vector3D;
        using CoeffVector = Vector<ScalarType, ElementType::DegreeOfFreedom>;
        using MDCellType = typename SolverType::MDCellType;
        using KSpaceFCGrid = typename SolverType::KSpaceFCGrid;
        using EigenValueGrid = GridStorage<Vector<ScalarType>>;
        constexpr static unsigned int Dim = Traits<MDCellType>::Dim;
        constexpr static unsigned int ElementVolume = 8;
    protected:
        SolverType solver;
        EigenValueGrid eigenvalues;
        CoeffVector diffCoeffX;
        CoeffVector diffCoeffY;
        CoeffVector diffCoeffZ;
    public:
        PhononDOS(MDCellType unitCell, Index3D superSize, const KSpaceFCGrid& forceConstants, Index3D gridDim);
        PhononDOS(const PhononDOS&) = default;
        PhononDOS(PhononDOS&&) noexcept = default;
        ~PhononDOS() = default;
        /* Operators */
        PhononDOS& operator=(PhononDOS obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType calcDOS(ScalarType freq) const;
        [[nodiscard]] ScalarType calcDOS(ScalarType freq, size_t band) const;
        [[nodiscard]] ScalarType calcDiffCapacityCv(ScalarType omegaW, ScalarType temperatureT);
        [[nodiscard]] ScalarType calcDiffHelmholtzF(ScalarType omegaW, ScalarType temperatureT);
        [[nodiscard]] ScalarType calcDiffEntropyS(ScalarType omegaW, ScalarType temperatureT);
        void swap(PhononDOS& __restrict obj) noexcept;
        /* Getter */
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return solver.getUnitCellDOF(); }
    protected:
        PhononDOS(MDCellType unitCell, Index3D superSize, Index3D gridDim);
        /* Operations */
        [[nodiscard]] ScalarType calcElemDOS(ScalarType freq, size_t band, Index3D index) const;
    };

    template<class ScalarType>
    PhononDOS<ScalarType>::PhononDOS(
            MDCellType unitCell, Index3D superSize, const KSpaceFCGrid& forceConstants, Index3D gridDim)
            : PhononDOS(std::move(unitCell), std::move(superSize), std::move(gridDim)) {
        eigenvalues.forIndexInGrid([this, &forceConstants](Index3D index) {
            const Index3D gridDim = eigenvalues.getDim();
            Vector3D qPoint{};
            for (unsigned int i = 0; i < Dim; ++i)
                qPoint[i] = ScalarType(index[i]) / ScalarType(gridDim[i]);
            auto fcMatrix = solver.interpolatePoint(qPoint, forceConstants);
            solver.toDynamicMatrix(fcMatrix);
            const auto eigen = SolverType::diagonalize(fcMatrix);
            eigenvalues(index) = solver.makeFreq(eigen);
        });
    }

    template<class ScalarType>
    PhononDOS<ScalarType>::PhononDOS(MDCellType unitCell, Index3D superSize, Index3D gridDim)
            : solver(std::move(unitCell), superSize), eigenvalues(gridDim) {
        for (unsigned i = 0; i < CoeffVector::SizeAtCompile; ++i) {
            diffCoeffX[i] = ElementType::dBase_dr(i);
            diffCoeffY[i] = ElementType::dBase_ds(i);
            diffCoeffZ[i] = ElementType::dBase_dt(i);
        }
    }

    template<class ScalarType>
    ScalarType PhononDOS<ScalarType>::calcDOS(ScalarType freq) const {
        ScalarType result = 0;
        for (size_t band = 0; band < solver.getNumBand(); ++band)
            result += calcDOS(freq, band);
        return result;
    }

    template<class ScalarType>
    ScalarType PhononDOS<ScalarType>::calcDOS(ScalarType freq, size_t band) const {
        ScalarType result = 0;
        eigenvalues.forIndexInGrid([this, freq, band, &result](Index3D index) {
            result += calcElemDOS(freq, band, index);
        });
        result /= ScalarType(eigenvalues.getSize() * ElementVolume);
        return result;
    }

    template<class ScalarType>
    ScalarType PhononDOS<ScalarType>::calcDiffCapacityCv(ScalarType omegaW, ScalarType temperatureT) {
        const ScalarType freq = omegaW * ScalarType(1 / (2 * M_PI));
        const ScalarType x = omegaW / temperatureT * 0.5;
        return calcDOS(freq) * square(x / sinh(x));
    }

    template<class ScalarType>
    ScalarType PhononDOS<ScalarType>::calcDiffHelmholtzF(ScalarType omegaW, ScalarType temperatureT) {
        const ScalarType freq = omegaW * ScalarType(1 / (2 * M_PI));
        const ScalarType dos = calcDOS(freq);
        const ScalarType zeroPointE = omegaW * 0.5;
        const ScalarType helmholtz1 = temperatureT * ln(ScalarType(1) - exp(-omegaW / temperatureT));
        return (zeroPointE + helmholtz1) * dos;
    }

    template<class ScalarType>
    ScalarType PhononDOS<ScalarType>::calcDiffEntropyS(ScalarType omegaW, ScalarType temperatureT) {
        const ScalarType freq = omegaW * ScalarType(1 / (2 * M_PI));
        const ScalarType dos = calcDOS(freq);
        const ScalarType x = omegaW / temperatureT;
        const ScalarType exp_x = exp(x);
        const ScalarType temp = exp_x - ScalarType(1);
        return (ln(temp) - x * exp_x / temp) * dos;
    }

    template<class ScalarType>
    void PhononDOS<ScalarType>::swap(PhononDOS& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        solver.swap(obj.solver);
        eigenvalues.swap(obj.eigenvalues);
        diffCoeffX.swap(obj.diffCoeffX);
        diffCoeffY.swap(obj.diffCoeffY);
        diffCoeffZ.swap(obj.diffCoeffZ);
    }

    template<class ScalarType>
    ScalarType PhononDOS<ScalarType>::calcElemDOS(ScalarType freq, size_t band, Index3D index) const {
        const Index3D gridDim = eigenvalues.getDim();
        Index3D index1{};
        for (unsigned int i = 0; i < Dim; ++i)
            index1[i] = (index[i] + 1) % gridDim[i];

        CoeffVector cornerFreq{};
        cornerFreq[0] = eigenvalues(index)[band];
        cornerFreq[1] = eigenvalues(index[0], index[1], index1[2])[band];
        cornerFreq[2] = eigenvalues(index[0], index1[1], index[2])[band];
        cornerFreq[3] = eigenvalues(index[0], index1[1], index1[2])[band];
        cornerFreq[4] = eigenvalues(index1[0], index[1], index[2])[band];
        cornerFreq[5] = eigenvalues(index1[0], index[1], index1[2])[band];
        cornerFreq[6] = eigenvalues(index1[0], index1[1], index[2])[band];
        cornerFreq[7] = eigenvalues(index1)[band];

        using Vector4D = Vector<ScalarType, 4>;
        const Vector4D plane{cornerFreq * diffCoeffX, cornerFreq * diffCoeffY, cornerFreq * diffCoeffZ, freq - mean(cornerFreq)};
        const auto head = plane.head(3);
        const auto cross = CubeCross<ScalarType>(plane);
        const ScalarType norm = head.norm();
        return cross.getArea() / norm;
    }
}
