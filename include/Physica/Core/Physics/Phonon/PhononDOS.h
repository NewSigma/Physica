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

#include "Physica/Core/Math/Calculus/PDE/FEM/Element/CuboidLinear.h"
#include "Physica/Core/Math/Geometry/CubeCross.h"
#include "PhononSolver.h"

namespace Physica {
    template<Scalar T>
    class PhononDOS {
    public:
        using SolverType = PhononSolver<T>;
        using ElementType = CuboidLinear<T>;
        using CoeffVector = DenseVector<T, ElementType::DegreeOfFreedom>;
        using MDCellType = SolverType::MDCellType;
        using KSpaceFCGrid = SolverType::KSpaceFCGrid;
        using EigenValueGrid = ArrayND<VectorND<T>, 3>;
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
        [[nodiscard]] T calcDOS(T freq) const;
        [[nodiscard]] T calcDOS(T freq, size_t band) const;
        [[nodiscard]] T calcDiffCapacityCv(T omegaW, T temperatureT);
        [[nodiscard]] T calcDiffHelmholtzF(T omegaW, T temperatureT);
        [[nodiscard]] T calcDiffEntropyS(T omegaW, T temperatureT);
        void swap(PhononDOS& __restrict obj) noexcept;
        /* Getter */
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return solver.getUnitCellDOF(); }
    protected:
        PhononDOS(MDCellType unitCell, Index3D superSize, Index3D gridDim);
        /* Operations */
        [[nodiscard]] T calcElemDOS(T freq, size_t band, Index3D index) const;
    };

    template<Scalar T>
    PhononDOS<T>::PhononDOS(
            MDCellType unitCell, Index3D superSize, const KSpaceFCGrid& forceConstants, Index3D gridDim)
            : PhononDOS(std::move(unitCell), std::move(superSize), std::move(gridDim)) {
        eigenvalues.forND([this, &forceConstants](VectorND<T>& eig, Index3D index) {
            const Index3D gridDim = eigenvalues.getShape();
            Vector3D<T> qPoint{};
            for (unsigned int i = 0; i < Dim; ++i)
                qPoint[i] = T(index[i]) / T(gridDim[i]);
            auto fcMatrix = solver.interpolatePoint(qPoint, forceConstants);
            solver.toDynamicMatrix(fcMatrix);
            const auto eigen = SolverType::diagonalize(fcMatrix);
            eig = solver.makeFreq(eigen);
        });
    }

    template<Scalar T>
    PhononDOS<T>::PhononDOS(MDCellType unitCell, Index3D superSize, Index3D gridDim)
            : solver(std::move(unitCell), superSize), eigenvalues(gridDim) {
        for (unsigned i = 0; i < CoeffVector::SizeAtCompile; ++i) {
            diffCoeffX[i] = ElementType::dBase_dr(i);
            diffCoeffY[i] = ElementType::dBase_ds(i);
            diffCoeffZ[i] = ElementType::dBase_dt(i);
        }
    }

    template<Scalar T>
    T PhononDOS<T>::calcDOS(T freq) const {
        T result = 0;
        for (size_t band = 0; band < solver.getNumBand(); ++band)
            result += calcDOS(freq, band);
        return result;
    }

    template<Scalar T>
    T PhononDOS<T>::calcDOS(T freq, size_t band) const {
        T result = 0;
        forND(eigenvalues.getShape(), [this, freq, band, &result](Index3D index) {
            result += calcElemDOS(freq, band, index);
        });
        result /= T(eigenvalues.getSize() * ElementVolume);
        return result;
    }

    template<Scalar T>
    T PhononDOS<T>::calcDiffCapacityCv(T omegaW, T temperatureT) {
        const T freq = omegaW * T(1 / (2 * M_PI));
        const T x = omegaW / temperatureT * 0.5;
        return calcDOS(freq) * square(x / sinh(x));
    }

    template<Scalar T>
    T PhononDOS<T>::calcDiffHelmholtzF(T omegaW, T temperatureT) {
        const T freq = omegaW * T(1 / (2 * M_PI));
        const T dos = calcDOS(freq);
        const T zeroPointE = omegaW * 0.5;
        const T helmholtz1 = temperatureT * ln(T(1) - exp(-omegaW / temperatureT));
        return (zeroPointE + helmholtz1) * dos;
    }

    template<Scalar T>
    T PhononDOS<T>::calcDiffEntropyS(T omegaW, T temperatureT) {
        const T freq = omegaW * T(1 / (2 * M_PI));
        const T dos = calcDOS(freq);
        const T x = omegaW / temperatureT;
        const T exp_x = exp(x);
        const T temp = exp_x - T(1);
        return (ln(temp) - x * exp_x / temp) * dos;
    }

    template<Scalar T>
    void PhononDOS<T>::swap(PhononDOS& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        solver.swap(obj.solver);
        eigenvalues.swap(obj.eigenvalues);
        diffCoeffX.swap(obj.diffCoeffX);
        diffCoeffY.swap(obj.diffCoeffY);
        diffCoeffZ.swap(obj.diffCoeffZ);
    }

    template<Scalar T>
    T PhononDOS<T>::calcElemDOS(T freq, size_t band, Index3D index) const {
        const Index3D shape = eigenvalues.getShape();
        Index3D index1{};
        for (unsigned int i = 0; i < Dim; ++i)
            index1[i] = (index[i] + 1) % shape[i];

        CoeffVector cornerFreq{};
        cornerFreq[0] = eigenvalues(index)[band];
        cornerFreq[1] = eigenvalues(index[0], index[1], index1[2])[band];
        cornerFreq[2] = eigenvalues(index[0], index1[1], index[2])[band];
        cornerFreq[3] = eigenvalues(index[0], index1[1], index1[2])[band];
        cornerFreq[4] = eigenvalues(index1[0], index[1], index[2])[band];
        cornerFreq[5] = eigenvalues(index1[0], index[1], index1[2])[band];
        cornerFreq[6] = eigenvalues(index1[0], index1[1], index[2])[band];
        cornerFreq[7] = eigenvalues(index1)[band];

        const Vector4D<T> plane{cornerFreq * diffCoeffX, cornerFreq * diffCoeffY, cornerFreq * diffCoeffZ, freq - cornerFreq.mean()};
        const auto head = plane.head(3);
        const auto cross = CubeCross<T>(plane);
        const T norm = head.norm();
        return cross.getArea() / norm;
    }
}
