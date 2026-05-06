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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "PhononSolverImpl/BHProjector.h"

namespace Physica {
    template<Scalar T>
    class PhononSolver {
        using This = PhononSolver<T>;
    public:
        using Tc = Complex<T>;
        using FFT3D = FFT<T, 3>;
        using MDCellType = MDCell<T>;
        using PositionMatrix = MDCellType::PositionMatrix;
        using RSpaceFCMat = FCProjector<T>::RSpaceFCMat;
        using RSpaceFCGrid = FCProjector<T>::RSpaceFCGrid;
        using KSpaceFCMat = DenseMatrix<Tc>;
        using KSpaceFCGrid = ArrayND<KSpaceFCMat, 3>;
        using EigenSolverType = EigenSolver<Tc>;
        using QPointGrid = ArrayND<EigenSolverType, 3>;
        using BornChargeArray = RSpaceEwald<T, true>::BornChargeArray;
        constexpr static unsigned int Dim = Traits<MDCellType>::Dim;
    private:
        using EwaldType = Ewald<T>;

        MDCellType unitCell;
        Index3D superSize;
    public:
        PhononSolver(MDCellType unitCell_, Index3D superSize_);
        PhononSolver(const This&) = default;
        PhononSolver(This&&) noexcept = default;
        ~PhononSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void removeTrans(KSpaceFCGrid& forceConstants, T translationPrec, size_t maxIteration) const;
        void projectTrans(RSpaceFCGrid& fcGrid, T maxAbsDot, size_t maxIteration) const;
        void projectTransRot(RSpaceFCGrid& fcGrid, T maxAbsDot, size_t maxIteration) const;
        void projectTransRotH(RSpaceFCGrid& fcGrid, T maxAbsDot, size_t maxIteration) const;

        [[nodiscard]] KSpaceFCGrid toKSpace(const RSpaceFCGrid& rSpaceGrid) const;
        [[nodiscard]] KSpaceFCMat interpolatePoint(Vector3D<T> qPoint, const KSpaceFCGrid& forceConstants) const;

        void applyNAC(RSpaceFCGrid& rSpaceGrid, const BornChargeArray& born) const;

        void toDynamicMatrix(KSpaceFCMat& forceConstant) const;
        void toDynamicMatrix(KSpaceFCGrid& forceConstants) const;
        [[nodiscard]] DenseMatrix<T> makeEigenVectors(const EigenSolverType& eigen) const;
        [[nodiscard]] DenseMatrix<T> makeEigenVectors(const QPointGrid& qPoints, Index3D qIndex) const;
        [[nodiscard]] MDCellType shiftAtom(const VectorND<T>& eigenVector, T distance);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const MDCellType& getUnitCell() const noexcept { return unitCell; }
        [[nodiscard]] size_t getNumUnitCellAtom() const noexcept { return unitCell.getNumParticle(); }
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return Dim * getNumUnitCellAtom(); }
        [[nodiscard]] size_t getNumBand() const noexcept { return getUnitCellDOF(); }
        [[nodiscard]] size_t getNumSuperCellAtom() const noexcept { return getNumUnitCellAtom() * getNumCell(); }
        [[nodiscard]] size_t getSuperCellDOF() const noexcept { return getUnitCellDOF() * getNumCell(); }
        [[nodiscard]] Index3D getSuperSize() const noexcept { return superSize; }
        [[nodiscard]] Index3D getForceConstantsGridSize() const noexcept { return FFT3D::rSizeToKSize(superSize); }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSize[0] * superSize[1] * superSize[2]; }
        /* Setters */
        void setUnitCell(MDCellType unitCell_) { unitCell = std::move(unitCell_); }
        /* Static members */
        [[nodiscard]] static EigenSolverType diagonalize(const KSpaceFCMat& dynamicMatrix);
        [[nodiscard]] static QPointGrid diagonalize(const KSpaceFCGrid& dynamicMatrixes);
        [[nodiscard]] static VectorND<T> makeFreq(const EigenSolverType& eigen);
        [[nodiscard]] static VectorND<T> makeFreq(const QPointGrid& qPoints, Index3D qIndex);
    private:
        T removeDriftForce(VectorND<T>& force) const;
    };

    template<Scalar T>
    PhononSolver<T>::PhononSolver(MDCellType unitCell_, Index3D superSize_)
            : unitCell(std::move(unitCell_))
            , superSize(superSize_) {}
    /**
     * Utilizing the iteration method to incorporate translational invariance while keep force constant matrix symmetric as introduced in [1].
     * 
     * References:
     * [1] Comput. Phys. Commun., 2009, 180(12), 2622-2633; https://doi.org/10.1016/j.cpc.2009.03.010
     */
    template<Scalar T>
    void PhononSolver<T>::removeTrans(
            KSpaceFCGrid& forceConstants,
            T translationPrec,
            size_t maxIteration) const {
        assert(translationPrec.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const size_t unitCellDOF = getUnitCellDOF();
        auto& fcMatrixes = forceConstants.asArray();
        FFT3D fft(superSize, PlanFlag::Estimate);
        auto rSpace = fft.getRSpace().flatten();
        auto kSpace = fft.getKSpace().flatten();

        VectorND<T> forceConst(getSuperCellDOF());
        KSpaceFCMat temp;
        T averageDrift = std::numeric_limits<T>::max();
        size_t iteration = 0;
        while (averageDrift > translationPrec) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Failed to apply phonon invariance in the given steps");
            iteration += 1;

            for (size_t row = 0; row < unitCellDOF; ++row) {
                for (size_t col = 0; col < unitCellDOF; ++col) {
                    for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                        auto& fcMatrix = fcMatrixes[i];
                        kSpace[i] = fcMatrix[row, col];
                    }
                    fft.invTransform();

                    const size_t shift = unitCellDOF;
                    for (size_t cell = 0; cell < getNumCell(); ++cell)
                        forceConst[col + cell * shift] = rSpace[cell];
                }
                averageDrift.toNextMean(row, removeDriftForce(forceConst));

                for (size_t col = 0; col < unitCellDOF; ++col) {
                    const size_t shift = unitCellDOF;
                    for (size_t cell = 0; cell < getNumCell(); ++cell)
                        rSpace[cell] = forceConst[col + cell * shift];
                    fft.transform();

                    for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                        auto& fcMatrix = fcMatrixes[i];
                        fcMatrix[row, col] = kSpace[i];
                    }
                }
            }

            for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                auto& fcMatrix = fcMatrixes[i];
                temp = (fcMatrix + fcMatrix.hermite()) * T(0.5);
                fcMatrix = temp;
            }
        }
    }

    template<Scalar T>
    void PhononSolver<T>::projectTrans(RSpaceFCGrid& fcGrid, T maxAbsDot, size_t maxIteration) const {
        assert(maxAbsDot.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const FCProjector<T> projector(superSize, getUnitCellDOF());
        auto fcVector = projector.toVector(fcGrid);
        T absDot = std::numeric_limits<T>::max();
        size_t iteration = 0;
        while (absDot > maxAbsDot) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Failed to apply phonon invariance in the given steps");
            iteration += 1;
            absDot = projector.projectSwap(fcVector);
            absDot = std::max(projector.projectTrans(fcVector), absDot);
        }
        projector.toGrid(fcVector, fcGrid);
    }

    template<Scalar T>
    void PhononSolver<T>::projectTransRot(RSpaceFCGrid& fcGrid, T maxAbsDot, size_t maxIteration) const {
        assert(maxAbsDot.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const BHProjector<T> projector(superSize, getUnitCellDOF(), unitCell);
        auto fcVector = projector.toVector(fcGrid);
        T absDot = std::numeric_limits<T>::max();
        size_t iteration = 0;
        while (absDot > maxAbsDot) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Failed to apply phonon invariance in the given steps");
            iteration += 1;
            absDot = projector.projectSwap(fcVector);
            absDot = std::max(projector.projectTrans(fcVector), absDot);
            absDot = std::max(projector.projectRot(fcVector), absDot);
        }
        projector.toGrid(fcVector, fcGrid);
    }

    template<Scalar T>
    void PhononSolver<T>::projectTransRotH(RSpaceFCGrid& fcGrid, T maxAbsDot, size_t maxIteration) const {
        assert(maxAbsDot.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const BHProjector<T> projector(superSize, getUnitCellDOF(), unitCell);
        auto fcVector = projector.toVector(fcGrid);
        T absDot = std::numeric_limits<T>::max();
        size_t iteration = 0;
        while (absDot > maxAbsDot) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Failed to apply phonon invariance in the given steps");
            iteration += 1;
            absDot = projector.projectSwap(fcVector);
            absDot = std::max(projector.projectTrans(fcVector), absDot);
            absDot = std::max(projector.projectRot(fcVector), absDot);
            absDot = std::max(projector.projectHuang(fcVector), absDot);
        }
        projector.toGrid(fcVector, fcGrid);
    }

    template<Scalar T>
    auto PhononSolver<T>::toKSpace(const RSpaceFCGrid& rSpaceGrid) const -> KSpaceFCGrid {
        assert(superSize == rSpaceGrid.getShape() && "[Error]: Super sizes do not match");
        assert((getUnitCellDOF() == rSpaceGrid[0, 0, 0].getRow()) && "[Error]: DOF do not match");
        const size_t unitCellDOF = getUnitCellDOF();
        KSpaceFCGrid kSpaceGrid(getForceConstantsGridSize(), unitCellDOF, unitCellDOF);
        FFT3D fft(superSize, PlanFlag::Estimate);
        for (size_t major = 0; major < unitCellDOF; ++major) {
            for (size_t minor = 0; minor < unitCellDOF; ++minor) {
                fft.getRSpace().forND([major, minor, &rSpaceGrid](T& elem, Index3D index) {
                    elem = rSpaceGrid[index].calcFromMajorMinor(major, minor);
                });
                fft.transform();

                fft.getKSpace().forND([major, minor, &kSpaceGrid](const Tc& elem, Index3D index) {
                    kSpaceGrid[index].refFromMajorMinor(major, minor) = elem;
                });
            }
        }
        return kSpaceGrid;
    }

    template<Scalar T>
    auto PhononSolver<T>::interpolatePoint(Vector3D<T> qPoint, const KSpaceFCGrid& forceConstants) const -> KSpaceFCMat {
        const Vector3D<T> qVector = unitCell.makeRepLattice().transpose() * qPoint;
        const size_t unitCellDOF = getUnitCellDOF();
        FFT3D fft(superSize, PlanFlag::Estimate);
        KSpaceFCMat result(unitCellDOF, unitCellDOF);
        for (size_t major = 0; major < unitCellDOF; ++major) {
            for (size_t minor = major; minor < unitCellDOF; ++minor) {
                auto& kSpace = fft.getKSpace();
                kSpace.forND([major, minor, &forceConstants](Tc& elem, Index3D index) {
                    elem = forceConstants[index].calcFromMajorMinor(major, minor);
                });
                fft.invTransform();

                Tc elem = 0;
                fft.getRSpace().forND([this, qVector, &elem](T& x, Index3D index) {
                    const auto& lattice = unitCell.getLattice();
                    T phase = 0;
                    T coeff = 1;
                    for (unsigned int i = 0; i < Dim; ++i) {
                        const ssize_t index_i = index[i];
                        const ssize_t dim_i = superSize[i];
                        const T factor = T(index_i > dim_i / 2 ? index_i - dim_i : (index_i));
                        const T phase_i = qVector * lattice.row(i) * factor;
                        const bool isOnWignerSeitzBoundary = (dim_i % 2 == 0) && (index_i == dim_i / 2);
                        if (isOnWignerSeitzBoundary)
                            coeff *= cos(phase_i);
                        else
                            phase += phase_i;
                    }
                    const auto factor = Tc::fromPhase(phase);
                    elem += x * coeff * factor;
                });
                result.refFromMajorMinor(major, minor) = elem;
                if (major != minor) [[likely]]
                    result.refFromMajorMinor(minor, major) = elem.conjugate();
            }
        }
        return result;
    }

    template<Scalar T>
    void PhononSolver<T>::applyNAC(RSpaceFCGrid& rSpaceGrid, const BornChargeArray& born) const {
        const T factor = reciprocal(unitCell.getVolume() * T(getNumCell() * PhyConst<AU>::vacuumDielectric));
        const auto repLatt = unitCell.makeRepLattice();
        rSpaceGrid.forND([this, factor, &repLatt, &born](T& rSpaceFC, Index3D index) {
            const bool isGammaPoint = index[0] == 0 && index[1] == 0 && index[2] == 0;
            Vector3D<T> qVector{};
            if (isGammaPoint)
                qVector = {1, 0, 0};
            else {
                Vector3D<T> qPoint{};
                for (unsigned int i = 0; i < Dim; ++i)
                    qPoint[i] = T(index[i]) / T(superSize[i]);
                qVector = repLatt.transpose() * qPoint;
                qVector.toUnit();
            }

            const size_t unitCellDOF = getUnitCellDOF();
            for (size_t major = 0; major < unitCellDOF; ++major) {
                const size_t atom1 = major / Dim;
                const size_t dir1 = major % Dim;
                const T projCharge1 = (born[atom1] * qVector).calc(dir1);
                for (size_t minor = major; minor < unitCellDOF; ++minor) {
                    const size_t atom2 = minor / Dim;
                    const size_t dir2 = minor % Dim;
                    const T projCharge2 = (born[atom2] * qVector).calc(dir2);
                    const T correction = factor * (projCharge1 * projCharge2);
                    rSpaceFC.refFromMajorMinor(major, minor) += correction;
                    if (minor != major)
                        rSpaceFC.refFromMajorMinor(minor, major) += correction;
                }
            }
        });
    }

    template<Scalar T>
    void PhononSolver<T>::toDynamicMatrix(KSpaceFCMat& forceConstant) const {
        const size_t unitCellDOF = getUnitCellDOF();
        for (size_t row = 0; row < unitCellDOF; ++row) {
            const size_t atom1 = row / Dim;
            const T mass1 = unitCell.getMass(atom1);
            for (size_t col = 0; col < unitCellDOF; ++col) {
                const size_t atom2 = col / Dim;
                const T mass2 = unitCell.getMass(atom2);
                const T repMass = reciprocal(sqrt(mass1 * mass2));
                forceConstant[row, col] *= repMass;
            }
        }
    }

    template<Scalar T>
    void PhononSolver<T>::toDynamicMatrix(KSpaceFCGrid& forceConstants) const {
        const size_t unitCellDOF = getUnitCellDOF();
        auto& fcMatrixes = forceConstants.flatten();
        for (size_t row = 0; row < unitCellDOF; ++row) {
            const size_t atom1 = row / Dim;
            const T mass1 = unitCell.getMass(atom1);
            for (size_t col = 0; col < unitCellDOF; ++col) {
                const size_t atom2 = col / Dim;
                const T mass2 = unitCell.getMass(atom2);
                const T repMass = reciprocal(sqrt(mass1 * mass2));

                for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                    auto& fcMatrix = fcMatrixes[i];
                    fcMatrix[row, col] *= repMass;
                }
            }
        }
    }

    template<Scalar T>
    DenseMatrix<T> PhononSolver<T>::makeEigenVectors(const EigenSolverType& eigen) const {
        DenseMatrix<T> result(getUnitCellDOF(), getUnitCellDOF());
        for (size_t i = 0; i < result.getCol(); ++i) {
            const auto fromCol = eigen.getRawEigenvectors().col(i);
            auto toCol = result.col(i);
            for (size_t j = 0; j < result.getRow(); ++j) {
                const T repSqrtMass = reciprocal(sqrt(unitCell.getMass(j / Dim)));
                toCol[j] = fromCol[j].real() * repSqrtMass;
            }
            toCol.toUnit(); // Optimize: repSqrtMass matrix is diagonal
        }
        return result;
    }

    template<Scalar T>
    DenseMatrix<T> PhononSolver<T>::makeEigenVectors(
            const QPointGrid& qPoints, Index3D qIndex) const {
        return makeEigenVectors(qPoints(qIndex));
    }

    template<Scalar T>
    auto PhononSolver<T>::shiftAtom(const VectorND<T>& eigenVector, T distance) -> MDCellType {
        assert(eigenVector.getLength() == getUnitCellDOF() && "[Error]: This is not a eigenVector of current cell");
        PositionMatrix shiftPos = unitCell.getPos() + eigenVector.reshape_row(getNumUnitCellAtom(), 3) * distance;
        return MDCellType(unitCell.getLattice(), std::move(shiftPos), unitCell.getMassVec());
    }

    template<Scalar T>
    void PhononSolver<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        unitCell.swap(obj.unitCell);
        superSize.swap(obj.superSize);
    }

    template<Scalar T>
    auto PhononSolver<T>::diagonalize(const KSpaceFCMat& dynamicMatrix) -> EigenSolverType {
        const size_t unitCellDOF = dynamicMatrix.getRow();
        auto eigen = EigenSolverType(unitCellDOF, true);
        eigen.compute(dynamicMatrix);
        eigen.sort();
        return eigen;
    }

    template<Scalar T>
    auto PhononSolver<T>::diagonalize(const KSpaceFCGrid& dynamicMatrixes) -> QPointGrid {
        const auto& matrixes = dynamicMatrixes.flatten();
        const size_t unitCellDOF = matrixes[0].getRow();
        QPointGrid qPoints(dynamicMatrixes.getShape(), unitCellDOF);
        for (size_t i = 0; i < matrixes.getLength(); ++i) {
            auto& eigen = qPoints.flatten()[i];
            eigen.compute(matrixes[i], true);
            eigen.sort();
        }
        return qPoints;
    }

    template<Scalar T>
    VectorND<T> PhononSolver<T>::makeFreq(const EigenSolverType& eigen) {
        const size_t unitCellDOF = eigen.getSize();
        VectorND<T> result(unitCellDOF);
        const VectorND<T> eigenvalues = eigen.getEigenvalues().reals();
        for (size_t i = 0; i < unitCellDOF; ++i)
            result[i] = eigenvalues[i].isNegative() ? -1 : 1;
        result = hadamard(result, sqrt(abs(eigenvalues)));
        result *= T(1 / (2 * M_PI));
        return result;
    }

    template<Scalar T>
    VectorND<T> PhononSolver<T>::makeFreq(const QPointGrid& qPoints, Index3D qIndex) {
        return makeFreq(qPoints(qIndex));
    }

    template<Scalar T>
    T PhononSolver<T>::removeDriftForce(VectorND<T>& force) const {
        const size_t superCellDOF = getSuperCellDOF();
        assert(force.getLength() == superCellDOF);
        const size_t numSuperAtom = superCellDOF / Dim;
        const T factor = reciprocal(T(numSuperAtom));
        T averageDrift = 0;
        for (size_t dim = 0; dim < Dim; ++dim) {
            T drift = 0;
            for (size_t i = dim; i < superCellDOF; i += Dim)
                drift += force[i];
            drift *= factor;
            for (size_t i = dim; i < superCellDOF; i += Dim)
                force[i] -= drift;
            averageDrift += abs(drift);
        }
        averageDrift /= T(Dim);
        return averageDrift;
    }
}
