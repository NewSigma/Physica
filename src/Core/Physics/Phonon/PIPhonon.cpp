/*
 * Copyright 2022-2024 WeiBo He.
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
#include "Physica/Core/Physics/Phonon/PIPhonon.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"

namespace Physica::Core {
    PIPhonon::PIPhonon(size_t numAtomUnitCell_,
                       size_t superSizeX_,
                       size_t superSizeY_,
                       size_t superSizeZ_)
            : fft({superSizeX_, superSizeY_, superSizeZ_}, PlanFlag::Estimate)
            , numAtomUnitCell(numAtomUnitCell_)
            , superSizeX(superSizeX_)
            , superSizeY(superSizeY_)
            , superSizeZ(superSizeZ_)
            , numSample(0) {
        force_corr.resize(getUnitCellDOF());
        momentum_corr.resize(getUnitCellDOF());
        const size_t numCell = getNumCell();
        for (size_t i = 0; i < force_corr.getLength(); ++i) {
            force_corr[i].resize(numCell, 0);
            momentum_corr[i].resize(numCell, 0);
        }

        kSpaceForceCorr.resize(numCell, getUnitCellDOF());
        kSpaceMomentumCorr.resize(numCell, getUnitCellDOF());

        normalModes.resize(numCell, getUnitCellDOF(), getUnitCellDOF());
    }

    PIPhonon& PIPhonon::operator=(PIPhonon obj) noexcept {
        swap(obj);
        return *this;
    }

    std::ostream& operator<<(std::ostream& os, const PIPhonon& phonon) {
        for (size_t r = 0; r < phonon.getUnitCellDOF(); ++r) {
            for (size_t c = r; c < phonon.getUnitCellDOF(); ++c) {
                const size_t offset = phonon.force_corr.toIndex1D(r, c);
                os << phonon.force_corr[offset] << phonon.momentum_corr[offset];
            }
        }
        return os;
    }

    std::istream& operator>>(std::istream& is, PIPhonon& phonon) {
        for (size_t r = 0; r < phonon.getUnitCellDOF(); ++r) {
            for (size_t c = r; c < phonon.getUnitCellDOF(); ++c) {
                const size_t offset = phonon.force_corr.toIndex1D(r, c);
                is >> phonon.force_corr[offset] >> phonon.momentum_corr[offset];
            }
        }
        return is;
    }

    void PIPhonon::compute() {
        toKSpace();
        applyTranslationInvariance(kSpaceForceCorr[0]);
        applyTranslationInvariance(kSpaceMomentumCorr[0]);

        DenseMatrix<ComplexType> buffer(getUnitCellDOF(), getUnitCellDOF());
        DenseMatrix<ComplexType> base(getUnitCellDOF(), getUnitCellDOF());
        EigenSolver<ComplexType> solver(buffer.getRow());
        for (size_t qPointId = 0; qPointId < kSpaceForceCorr.getLength(); ++qPointId) {
            const bool isGammaPoint = qPointId == 0;
            if (isGammaPoint) {
                const DenseSymmMatrix<ScalarType> m = toRealMatrix(kSpaceMomentumCorr[qPointId]);
                const Schur<ScalarType> schur(m, true);
                for (size_t i = 0; i < base.getColumn(); ++i) {
                    const ScalarType eigenvalue = schur.getMatrixT().diag().calc(i);
                    if (eigenvalue > ScalarType(ConsiderAsZeroThrehold))
                        base.col(i) = schur.getMatrixU().col(i).asVector() * reciprocal(sqrt(eigenvalue));
                    else
                        base.col(i).asVector() = ScalarType(0);
                }
            }
            else {
                const Schur<ComplexType> schur(kSpaceMomentumCorr[qPointId], true);
                for (size_t i = 0; i < base.getColumn(); ++i) {
                    const ScalarType eigenvalue = schur.getMatrixT().diag().calc(i).getReal();
                    if (eigenvalue > ScalarType(ConsiderAsZeroThrehold))
                        base.col(i) = schur.getMatrixU().col(i).asVector() * reciprocal(sqrt(eigenvalue));
                    else
                        base.col(i).asVector() = ScalarType(0);
                }
            }
            buffer = base.hermite() * kSpaceForceCorr[qPointId].asMatrix();
            buffer *= base;
            
            solver.compute(buffer, true);
            solver.sort();
            normalModes[qPointId] = base * solver.getEigenvectors();
            for (size_t i = 0; i < buffer.getColumn(); ++i)
                kSpaceMomentumCorr[qPointId](i, i) = solver.getEigenvalues()[i];
        }
    }

    void PIPhonon::swap(PIPhonon& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        fft.swap(obj.fft);
        std::swap(numAtomUnitCell, obj.numAtomUnitCell);
        std::swap(superSizeX, obj.superSizeX);
        std::swap(superSizeY, obj.superSizeY);
        std::swap(superSizeZ, obj.superSizeZ);
        force_corr.swap(obj.force_corr);
        momentum_corr.swap(obj.momentum_corr);
        std::swap(numSample, obj.numSample);
        kSpaceForceCorr.swap(obj.kSpaceForceCorr);
        kSpaceMomentumCorr.swap(obj.kSpaceMomentumCorr);
        normalModes.swap(obj.normalModes);
    }

    void PIPhonon::toKSpace() {
        const size_t dof = getUnitCellDOF();
        for (size_t r = 0; r < dof; ++r) {
            for (size_t c = r; c < dof; ++c) {
                const size_t offset_corr = force_corr.toIndex1D(r, c);
                fft.getRSpace().flatten() = force_corr[offset_corr];
                fft.transform();
                for (size_t cell = 0; cell < getNumCell(); ++cell)
                    kSpaceForceCorr[cell](r, c) = fft.getKSpace().flatten()[cell];

                fft.getRSpace().flatten() = momentum_corr[offset_corr];
                fft.transform();
                for (size_t cell = 0; cell < getNumCell(); ++cell)
                    kSpaceMomentumCorr[cell](r, c) = fft.getKSpace().flatten()[cell];
            }
        }
    }

    void PIPhonon::applyTranslationInvariance(DenseHermiteMatrix<ComplexType>& target) {
        [[maybe_unused]] const bool isRealMatrix = &target == &kSpaceForceCorr[0] || &target == &kSpaceMomentumCorr[0];
        assert(isRealMatrix);
        DenseSymmMatrix<ScalarType> buffer(target.getRow(), target.getColumn());
        for (size_t c = 0; c < target.getColumn(); ++c) {
            for (size_t r = 0; r <= c; ++r) {
                ScalarType term1 = 0;
                for (size_t i = 0; i < getNumAtomUnitCell(); ++i)
                    for (size_t j = 0; j < getNumAtomUnitCell(); ++j)
                        term1 += target.calc(3 * i + r % 3U, 3 * j + c % 3U).getReal();
                term1 /= ScalarType(getNumAtomUnitCell() * getNumAtomUnitCell());

                ScalarType term2 = 0;
                for (size_t i = 0; i < getNumAtomUnitCell(); ++i)
                    term2 += target.calc(r, 3 * i + c % 3U).getReal() + target.calc(3 * i + r % 3U, c).getReal();
                term2 /= ScalarType(getNumAtomUnitCell());
                buffer(r, c) = term1 - term2;
            }
        }

        for (size_t c = 0; c < target.getColumn(); ++c)
            for (size_t r = 0; r <= c; ++r)
                target(r, c) += buffer(r, c);
    }
}
