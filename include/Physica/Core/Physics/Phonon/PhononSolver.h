/*
 * Copyright 2023-2024 WeiBo He.
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

#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "PhononSolverImpl/BHProjector.h"

namespace Physica::Core {
    template<class ScalarType>
    class PhononSolver {
        using This = PhononSolver<ScalarType>;
    public:
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using FFT3D = FFT<ScalarType, 3>;
        using MDCellType = MDCell<ScalarType>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using Index3D = typename FCProjector<ScalarType>::Index3D;
        using RSpaceFCMat = typename FCProjector<ScalarType>::RSpaceFCMat;
        using RSpaceFCGrid = typename FCProjector<ScalarType>::RSpaceFCGrid;
        using KSpaceFCMat = DenseMatrix<ComplexType>;
        using KSpaceFCGrid = GridStorage<KSpaceFCMat>;
        using EigenSolverType = EigenSolver<ComplexType>;
        using QPointGrid = GridStorage<EigenSolverType>;
        using BornChargeArray = typename RSpaceEwald<ScalarType, true>::BornChargeArray;
        constexpr static unsigned int Dim = Traits<MDCellType>::Dim;
    private:
        using EwaldType = Ewald<ScalarType>;

        MDCellType unitCell;
        Index3D superSize;
    public:
        PhononSolver(MDCellType unitCell_, Index3D superSize_);
        PhononSolver(const PhononSolver&) = default;
        PhononSolver(PhononSolver&&) noexcept = default;
        ~PhononSolver() = default;
        /* Operators */
        PhononSolver& operator=(PhononSolver obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void removeTrans(KSpaceFCGrid& forceConstants, ScalarType translationPrec, size_t maxIteration) const;
        void projectTrans(RSpaceFCGrid& fcGrid, ScalarType maxAbsDot, size_t maxIteration) const;
        void projectTransRot(RSpaceFCGrid& fcGrid, ScalarType maxAbsDot, size_t maxIteration) const;
        void projectTransRotH(RSpaceFCGrid& fcGrid, ScalarType maxAbsDot, size_t maxIteration) const;

        [[nodiscard]] KSpaceFCGrid toKSpace(const RSpaceFCGrid& rSpaceGrid) const;
        [[nodiscard]] KSpaceFCMat interpolatePoint(Vector3D qPoint, const KSpaceFCGrid& forceConstants) const;

        void applyNAC(RSpaceFCGrid& rSpaceGrid, const BornChargeArray& born) const;

        void toDynamicMatrix(KSpaceFCMat& forceConstant) const;
        void toDynamicMatrix(KSpaceFCGrid& forceConstants) const;
        [[nodiscard]] DenseMatrix<ScalarType> makeEigenVectors(const EigenSolverType& eigen) const;
        [[nodiscard]] inline DenseMatrix<ScalarType> makeEigenVectors(const QPointGrid& qPoints, Index3D qIndex) const;
        [[nodiscard]] MDCellType shiftAtom(const Vector<ScalarType>& eigenVector, ScalarType distance);
        void swap(PhononSolver& __restrict obj) noexcept;
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
        [[nodiscard]] static Vector<ScalarType> makeFreq(const EigenSolverType& eigen);
        [[nodiscard]] static inline Vector<ScalarType> makeFreq(const QPointGrid& qPoints, Index3D qIndex);
    private:
        ScalarType removeDriftForce(Vector<ScalarType>& force) const;
    };

    template<class ScalarType>
    PhononSolver<ScalarType>::PhononSolver(MDCellType unitCell_, Index3D superSize_)
            : unitCell(std::move(unitCell_))
            , superSize(superSize_) {}
    /**
     * Utilizing the iteration method to incorporate translational invariance while keep force constant matrix symmetric as introduced in [1].
     * 
     * References:
     * [1] Comput. Phys. Commun., 2009, 180(12), 2622-2633; https://doi.org/10.1016/j.cpc.2009.03.010
     */
    template<class ScalarType>
    void PhononSolver<ScalarType>::removeTrans(
            KSpaceFCGrid& forceConstants,
            ScalarType translationPrec,
            size_t maxIteration) const {
        assert(translationPrec.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const size_t unitCellDOF = getUnitCellDOF();
        auto& fcMatrixes = forceConstants.asArray();
        FFT3D fft(superSize, PlanFlag::Estimate);
        auto rSpace = fft.getRSpace().flatten();
        auto kSpace = fft.getKSpace().flatten();

        Vector<ScalarType> forceConst(getSuperCellDOF());
        KSpaceFCMat temp;
        ScalarType averageDrift = std::numeric_limits<ScalarType>::max();
        size_t iteration = 0;
        while (averageDrift > translationPrec) {
            if (iteration == maxIteration) [[unlikely]]
                throw BadConvergenceException("[Error]: Failed to apply phonon invariance in the given steps");
            iteration += 1;

            for (size_t row = 0; row < unitCellDOF; ++row) {
                for (size_t col = 0; col < unitCellDOF; ++col) {
                    for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                        auto& fcMatrix = fcMatrixes[i];
                        kSpace[i] = fcMatrix(row, col);
                    }
                    fft.invTransform();

                    const size_t shift = unitCellDOF;
                    for (size_t cell = 0; cell < getNumCell(); ++cell)
                        forceConst[col + cell * shift] = rSpace[cell];
                }
                toNextMean(averageDrift, row, removeDriftForce(forceConst));

                for (size_t col = 0; col < unitCellDOF; ++col) {
                    const size_t shift = unitCellDOF;
                    for (size_t cell = 0; cell < getNumCell(); ++cell)
                        rSpace[cell] = forceConst[col + cell * shift];
                    fft.transform();

                    for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                        auto& fcMatrix = fcMatrixes[i];
                        fcMatrix(row, col) = kSpace[i];
                    }
                }
            }

            for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                auto& fcMatrix = fcMatrixes[i];
                temp = (fcMatrix + fcMatrix.hermite()) * ScalarType(0.5);
                fcMatrix = temp;
            }
        }
    }

    template<class ScalarType>
    void PhononSolver<ScalarType>::projectTrans(RSpaceFCGrid& fcGrid, ScalarType maxAbsDot, size_t maxIteration) const {
        assert(maxAbsDot.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const FCProjector<ScalarType> projector(superSize, getUnitCellDOF());
        auto fcVector = projector.toVector(fcGrid);
        ScalarType absDot = std::numeric_limits<ScalarType>::max();
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

    template<class ScalarType>
    void PhononSolver<ScalarType>::projectTransRot(RSpaceFCGrid& fcGrid, ScalarType maxAbsDot, size_t maxIteration) const {
        assert(maxAbsDot.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const BHProjector<ScalarType> projector(superSize, getUnitCellDOF(), unitCell);
        auto fcVector = projector.toVector(fcGrid);
        ScalarType absDot = std::numeric_limits<ScalarType>::max();
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

    template<class ScalarType>
    void PhononSolver<ScalarType>::projectTransRotH(RSpaceFCGrid& fcGrid, ScalarType maxAbsDot, size_t maxIteration) const {
        assert(maxAbsDot.isPositive() && "[Error]: Invalid param");
        assert(maxIteration > 0 && "[Error]: Zero max iteration does nothing");
        const BHProjector<ScalarType> projector(superSize, getUnitCellDOF(), unitCell);
        auto fcVector = projector.toVector(fcGrid);
        ScalarType absDot = std::numeric_limits<ScalarType>::max();
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

    template<class ScalarType>
    typename PhononSolver<ScalarType>::KSpaceFCGrid PhononSolver<ScalarType>::toKSpace(const RSpaceFCGrid& rSpaceGrid) const {
        assert(superSize == rSpaceGrid.getDim() && "[Error]: Super sizes do not match");
        assert(getUnitCellDOF() == rSpaceGrid(0, 0, 0).getRow() && "[Error]: DOF do not match");
        const size_t unitCellDOF = getUnitCellDOF();
        KSpaceFCGrid kSpaceGrid(getForceConstantsGridSize(), unitCellDOF, unitCellDOF);
        FFT3D fft(superSize, PlanFlag::Estimate);
        for (size_t major = 0; major < unitCellDOF; ++major) {
            for (size_t minor = 0; minor < unitCellDOF; ++minor) {
                auto& rSpaceFFT = fft.getRSpace();
                rSpaceFFT.forIndexInGrid([major, minor, &rSpaceGrid, &rSpaceFFT](Index3D index) {
                    rSpaceFFT(index) = rSpaceGrid(index).calcFromMajorMinor(major, minor);
                });
                fft.transform();

                const auto& kSpaceFFT = fft.getKSpace();
                kSpaceFFT.forIndexInGrid([major, minor, &kSpaceGrid, &kSpaceFFT](Index3D index) {
                    kSpaceGrid(index).refFromMajorMinor(major, minor) = kSpaceFFT(index);
                });
            }
        }
        return kSpaceGrid;
    }

    template<class ScalarType>
    typename PhononSolver<ScalarType>::KSpaceFCMat PhononSolver<ScalarType>::interpolatePoint(
            Vector3D qPoint, const KSpaceFCGrid& forceConstants) const {
        const Vector3D qVector = unitCell.makeRepLattice().transpose() * qPoint;
        const size_t unitCellDOF = getUnitCellDOF();
        FFT3D fft(superSize, PlanFlag::Estimate);
        KSpaceFCMat result(unitCellDOF, unitCellDOF);
        for (size_t major = 0; major < unitCellDOF; ++major) {
            for (size_t minor = major; minor < unitCellDOF; ++minor) {
                auto& kSpace = fft.getKSpace();
                kSpace.forIndexInGrid([major, minor, &kSpace, &forceConstants](Index3D index) {
                    kSpace(index) = forceConstants(index).calcFromMajorMinor(major, minor);
                });
                fft.invTransform();

                ComplexType elem = 0;
                auto& rSpace = fft.getRSpace();
                rSpace.forIndexInGrid([this, qVector, &rSpace, &elem](Index3D index) {
                    const auto& lattice = unitCell.getLattice();
                    const Index3D rSpaceDim = rSpace.getDim();
                    ScalarType phase = 0;
                    ScalarType coeff = 1;
                    for (unsigned int i = 0; i < Dim; ++i) {
                        const ssize_t index_i = index[i];
                        const ssize_t dim_i = rSpaceDim[i];
                        const ScalarType factor = ScalarType(index_i > dim_i / 2 ? index_i - dim_i : (index_i));
                        const ScalarType phase_i = qVector * lattice.row(i).asVector() * factor;
                        const bool isOnWignerSeitzBoundary = (dim_i % 2 == 0) && (index_i == dim_i / 2);
                        if (isOnWignerSeitzBoundary)
                            coeff *= cos(phase_i);
                        else
                            phase += phase_i;
                    }
                    const auto factor = ComplexType::fromPhase(phase);
                    elem += rSpace(index) * coeff * factor;
                });
                result.refFromMajorMinor(major, minor) = elem;
                if (major != minor) [[likely]]
                    result.refFromMajorMinor(minor, major) = elem.conjugate();
            }
        }
        return result;
    }

    template<class ScalarType>
    void PhononSolver<ScalarType>::applyNAC(RSpaceFCGrid& rSpaceGrid, const BornChargeArray& born) const {
        const ScalarType factor = reciprocal(unitCell.getVolume() * ScalarType(getNumCell() * PhyConst<AU>::vacuumDielectric));
        const auto repLatt = unitCell.makeRepLattice();
        rSpaceGrid.forIndexInGrid([this, factor, &repLatt, &rSpaceGrid, &born](Index3D index) {
            const bool isGammaPoint = index[0] == 0 && index[1] == 0 && index[2] == 0;
            Vector3D qVector{};
            if (isGammaPoint)
                qVector = {1, 0, 0};
            else {
                Vector3D qPoint{};
                for (unsigned int i = 0; i < Dim; ++i)
                    qPoint[i] = ScalarType(index[i]) / ScalarType(superSize[i]);
                qVector = repLatt.transpose() * qPoint;
                qVector.toUnit();
            }

            auto& rSpaceFC = rSpaceGrid(index);
            const size_t unitCellDOF = getUnitCellDOF();
            for (size_t major = 0; major < unitCellDOF; ++major) {
                const size_t atom1 = major / Dim;
                const size_t dir1 = major % Dim;
                const ScalarType projCharge1 = (born[atom1] * qVector).calc(dir1);
                for (size_t minor = major; minor < unitCellDOF; ++minor) {
                    const size_t atom2 = minor / Dim;
                    const size_t dir2 = minor % Dim;
                    const ScalarType projCharge2 = (born[atom2] * qVector).calc(dir2);
                    const ScalarType correction = factor * (projCharge1 * projCharge2);
                    rSpaceFC.refFromMajorMinor(major, minor) += correction;
                    if (minor != major)
                        rSpaceFC.refFromMajorMinor(minor, major) += correction;
                }
            }
        });
    }

    template<class ScalarType>
    void PhononSolver<ScalarType>::toDynamicMatrix(KSpaceFCMat& forceConstant) const {
        const size_t unitCellDOF = getUnitCellDOF();
        for (size_t row = 0; row < unitCellDOF; ++row) {
            const size_t atom1 = row / Dim;
            const ScalarType mass1 = unitCell.getMass(atom1);
            for (size_t col = 0; col < unitCellDOF; ++col) {
                const size_t atom2 = col / Dim;
                const ScalarType mass2 = unitCell.getMass(atom2);
                const ScalarType repMass = reciprocal(sqrt(mass1 * mass2));
                forceConstant(row, col) *= repMass;
            }
        }
    }

    template<class ScalarType>
    void PhononSolver<ScalarType>::toDynamicMatrix(KSpaceFCGrid& forceConstants) const {
        const size_t unitCellDOF = getUnitCellDOF();
        auto& fcMatrixes = forceConstants.flatten();
        for (size_t row = 0; row < unitCellDOF; ++row) {
            const size_t atom1 = row / Dim;
            const ScalarType mass1 = unitCell.getMass(atom1);
            for (size_t col = 0; col < unitCellDOF; ++col) {
                const size_t atom2 = col / Dim;
                const ScalarType mass2 = unitCell.getMass(atom2);
                const ScalarType repMass = reciprocal(sqrt(mass1 * mass2));

                for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                    auto& fcMatrix = fcMatrixes[i];
                    fcMatrix(row, col) *= repMass;
                }
            }
        }
    }

    template<class ScalarType>
    DenseMatrix<ScalarType> PhononSolver<ScalarType>::makeEigenVectors(const EigenSolverType& eigen) const {
        DenseMatrix<ScalarType> result(getUnitCellDOF(), getUnitCellDOF());
        for (size_t i = 0; i < result.getColumn(); ++i) {
            const auto fromCol = eigen.getRawEigenvectors().col(i);
            auto toCol = result.col(i);
            for (size_t j = 0; j < result.getRow(); ++j) {
                const ScalarType repSqrtMass = reciprocal(sqrt(unitCell.getMass(j / Dim)));
                toCol[j] = fromCol[j].getReal() * repSqrtMass;
            }
            toCol.toUnit(); // Optimize: repSqrtMass matrix is diagonal
        }
        return result;
    }

    template<class ScalarType>
    inline DenseMatrix<ScalarType> PhononSolver<ScalarType>::makeEigenVectors(
            const QPointGrid& qPoints, Index3D qIndex) const {
        return makeEigenVectors(qPoints(qIndex));
    }

    template<class ScalarType>
    typename PhononSolver<ScalarType>::MDCellType
    PhononSolver<ScalarType>::shiftAtom(const Vector<ScalarType>& eigenVector, ScalarType distance) {
        assert(eigenVector.getLength() == getUnitCellDOF() && "[Error]: This is not a eigenVector of current cell");
        PositionMatrix shiftPos = unitCell.getPos() + eigenVector.reshape_row(getNumUnitCellAtom(), 3) * distance;
        return MDCellType(unitCell.getLattice(), std::move(shiftPos), unitCell.getMassVec());
    }

    template<class ScalarType>
    void PhononSolver<ScalarType>::swap(PhononSolver& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        unitCell.swap(obj.unitCell);
        superSize.swap(obj.superSize);
    }

    template<class ScalarType>
    typename PhononSolver<ScalarType>::EigenSolverType
    PhononSolver<ScalarType>::diagonalize(const KSpaceFCMat& dynamicMatrix) {
        const size_t unitCellDOF = dynamicMatrix.getRow();
        auto eigen = EigenSolverType(unitCellDOF);
        eigen.compute(dynamicMatrix, true);
        eigen.sort();
        return eigen;
    }

    template<class ScalarType>
    typename PhononSolver<ScalarType>::QPointGrid PhononSolver<ScalarType>::diagonalize(
            const KSpaceFCGrid& dynamicMatrixes) {
        const auto& matrixes = dynamicMatrixes.flatten();
        const size_t unitCellDOF = matrixes[0].getRow();
        QPointGrid qPoints(dynamicMatrixes.getDim(), unitCellDOF);
        for (size_t i = 0; i < matrixes.getLength(); ++i) {
            auto& eigen = qPoints.flatten()[i];
            eigen.compute(matrixes[i], true);
            eigen.sort();
        }
        return qPoints;
    }

    template<class ScalarType>
    Vector<ScalarType> PhononSolver<ScalarType>::makeFreq(const EigenSolverType& eigen) {
        const size_t unitCellDOF = eigen.getSize();
        Vector<ScalarType> result(unitCellDOF);
        const Vector<ScalarType> eigenvalues = toRealVector(eigen.getEigenvalues());
        for (size_t i = 0; i < unitCellDOF; ++i)
            result[i] = eigenvalues[i].isNegative() ? -1 : 1;
        result = hadamard(result, sqrt(abs(eigenvalues)));
        result *= ScalarType(1 / (2 * M_PI));
        return result;
    }

    template<class ScalarType>
    inline Vector<ScalarType> PhononSolver<ScalarType>::makeFreq(const QPointGrid& qPoints, Index3D qIndex) {
        return makeFreq(qPoints(qIndex));
    }

    template<class ScalarType>
    ScalarType PhononSolver<ScalarType>::removeDriftForce(Vector<ScalarType>& force) const {
        const size_t superCellDOF = getSuperCellDOF();
        assert(force.getLength() == superCellDOF);
        const size_t numSuperAtom = superCellDOF / Dim;
        const ScalarType factor = reciprocal(ScalarType(numSuperAtom));
        ScalarType averageDrift = 0;
        for (size_t dim = 0; dim < Dim; ++dim) {
            ScalarType drift = 0;
            for (size_t i = dim; i < superCellDOF; i += Dim)
                drift += force[i];
            drift *= factor;
            for (size_t i = dim; i < superCellDOF; i += Dim)
                force[i] -= drift;
            averageDrift += abs(drift);
        }
        averageDrift /= ScalarType(Dim);
        return averageDrift;
    }
}
