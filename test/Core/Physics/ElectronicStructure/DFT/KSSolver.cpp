/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/XCProvider/LDA.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using ComplexType = Complex<ScalarType>;
using RandomSource = Random<>;
using LDAProvider = LDA<ScalarType, LDAType::HL, false>;
using SolverType = KSSolver<ScalarType, LDAProvider>;
using BasisType = PlainWaveBasis<ScalarType>;
using DensityType = DensityGrid<ScalarType, false>;
using MixerType = ChargeMixer<ScalarType, false>;
using HamiltonType = KSHamilton<ScalarType, false>;
using BandType = BandGrid<ScalarType, false>;
using CrystalCellType = CrystalCell<ScalarType>;
using LatticeMatrix = CrystalCellType::LatticeMatrix;
using PositionMatrix = CrystalCellType::PositionMatrix;

template class Physica::KSSolver<ScalarType, LDAProvider>;
template void Physica::KSSolver<ScalarType, LDAProvider>::initWaveFunc<RandomSource>(const SpinPair<BasisType, false>&);

namespace {
    constexpr double Precision = 1E-12;

    CrystalCellType makeUnitCell() {
        const LatticeMatrix lattice{4, 0, 0, 0, 4, 0, 0, 0, 4};
        const PositionMatrix pos{0.5, 0.5, 0.5};
        return CrystalCellType({lattice, pos, CrystalCellType::Type::Direct}, {1});
    }

    void testForKInBasis(const LatticeMatrix& repLatt) {
        const Index3D dim{3, 4, 5};
        const size_t numK = dim[0] * dim[1] * dim[2];
        Array<Vector3D<ScalarType>> kAlone(numK);
        Array<Vector3D<ScalarType>> kWithIndex(numK);

        size_t count = 0;
        BasisType::forKInBasis(repLatt, dim, [&kAlone, &count](Vector3D<ScalarType> waveK) {
            kAlone[count] = waveK;
            ++count;
        });
        expect(count == numK);

        count = 0;
        bool isIndexConsistent = true;
        BasisType::forKIndexInBasis(repLatt, dim, [&](Vector3D<ScalarType> waveK, Index3D index) {
            kWithIndex[count] = waveK;
            isIndexConsistent = isIndexConsistent && waveK == BasisType::makeWaveVector(repLatt, index, dim);
            ++count;
        });
        expect(count == numK);
        expect(isIndexConsistent);

        bool isSameGrid = true;
        for (size_t i = 0; i < numK; ++i)
            isSameGrid = isSameGrid && kAlone[i] == kWithIndex[i];
        expect(isSameGrid);
    }

    void testDensityInterpolation() {
        const Index3D dim{4, 5, 6};
        DensityType density(dim);
        auto& grid = density.getTotalDensity();
        const auto toFraction = [&dim](size_t index, int axis) {
            return ScalarType(index) / ScalarType(dim[axis]);
        };
        // Multilinear field: trilinear interpolation reproduces it exactly
        forND(dim, [&](Index3D index) {
            const ScalarType x = toFraction(index[0], 0);
            const ScalarType y = toFraction(index[1], 1);
            const ScalarType z = toFraction(index[2], 2);
            grid[index] = ScalarType(1) + ScalarType(2) * x - ScalarType(3) * y + ScalarType(5) * x * z;
        });
        bool isReproduced = true;
        for (size_t i = 0; i + 1 < dim[0]; ++i) {
            for (size_t j = 0; j + 1 < dim[1]; ++j) {
                for (size_t k = 0; k + 1 < dim[2]; ++k) {
                    const ScalarType x = ScalarType(i + 0.37) / ScalarType(dim[0]);
                    const ScalarType y = ScalarType(j + 0.61) / ScalarType(dim[1]);
                    const ScalarType z = ScalarType(k + 0.19) / ScalarType(dim[2]);
                    const ScalarType answer = ScalarType(1) + ScalarType(2) * x - ScalarType(3) * y + ScalarType(5) * x * z;
                    isReproduced = isReproduced && scalarNear(density(SpinState::Up, {x, y, z}), answer, Precision);
                }
            }
        }
        expect(isReproduced);
        // Constant field: the wrap around of the last cell must hit the first one
        grid = ScalarType(7);
        bool isWrapped = true;
        for (size_t i = 0; i < dim[0]; ++i) {
            const ScalarType x = ScalarType(i + 0.9) / ScalarType(dim[0]);
            const ScalarType y = ScalarType(0.5) / ScalarType(dim[1]);
            const ScalarType z = ScalarType(0.5) / ScalarType(dim[2]);
            isWrapped = isWrapped && scalarNear(density(SpinState::Up, {x, y, z}), ScalarType(7), Precision);
        }
        expect(isWrapped);
    }

    void testBandGridSwap(const LatticeMatrix& repLatt) {
        BandType band1(repLatt, 2, 2, 2, 4, 4);
        BandType band2(repLatt, 3, 1, 1, 2, 3);
        band1.swap(band2);

        const Index3D shape1 = band1.getKPointGrid().getShape();
        const Index3D shape2 = band2.getKPointGrid().getShape();
        expect(shape1[0] == 3 && shape1[1] == 1 && shape1[2] == 1);
        expect(shape2[0] == 2 && shape2[1] == 2 && shape2[2] == 2);
        expect(band1.getNumBand() == 3);
        expect(band2.getNumBand() == 4);
    }

    void testKSHamilton(const CrystalCellType& cell) {
        constexpr ScalarType cutEnergyPsi = 1;
        constexpr ScalarType cutEnergyRho = 5;
        constexpr double kineticFactor = PhyConst<AU>::reducedPlanck * PhyConst<AU>::reducedPlanck
                                       / PhyConst<AU>::electronMass * 0.5;
        const Vector3D<ScalarType> kPoint{0, 0, 0};
        HamiltonType hamilton(cell, cutEnergyPsi, cutEnergyRho);
        const size_t numBasis = hamilton.getNumBasis();

        hamilton.makeFreeHamilton(kPoint);
        Array<ScalarType> freeDiagonal(numBasis);
        for (size_t i = 0; i < numBasis; ++i)
            freeDiagonal[i] = hamilton[SpinState::Up].calc(i, i).real();

        bool isKinetic = true;
        size_t count = 0;
        const Index3D basisDim = BasisType::makeGridDim(cutEnergyPsi, hamilton.getRepLattice());
        BasisType::forKInBasis(hamilton.getRepLattice(), basisDim, [&](Vector3D<ScalarType> waveK) {
            const ScalarType kineticE = (kPoint + waveK).squaredNorm() * ScalarType(kineticFactor);
            isKinetic = isKinetic && scalarNear(freeDiagonal[count], kineticE, Precision);
            ++count;
        });
        expect(count == numBasis);
        expect(isKinetic);

        /* The vanishing G = 0 component of the ion potential leaves the diagonal untouched */
        hamilton.makeNearFreeHamilton(kPoint);
        bool isUntouched = true;
        for (size_t i = 0; i < numBasis; ++i) {
            const ScalarType diagonal = hamilton[SpinState::Up].calc(i, i).real();
            isUntouched = isUntouched && scalarNear(diagonal, freeDiagonal[i], Precision);
        }
        expect(isUntouched);
        expect(isUntouched);
    }

    void testChargeMixer(const LatticeMatrix& repLatt) {
        const Index3D dim{4, 4, 4};
        MixerType mixer(repLatt, dim);
        DensityType input(dim);
        DensityType output(dim);
        DensityType mixed(dim);
        input.getTotalDensity().flatten() = ScalarType(2);
        output.getTotalDensity().flatten() = ScalarType(5);

        mixer.fetchInputDensity(input);
        mixer.mix(output);
        expect(scalarNear(mixer.getResidule().getTotalDensity().flatten()[0], ScalarType(3), Precision));

        /* A constant residual only lives at the Gamma point, where the mixing kernel vanishes */
        mixer.mix(size_t(0), mixed);
        const auto& result = mixed.getTotalDensity().flatten();
        bool isMixed = true;
        for (size_t i = 0; i < result.getLength(); ++i)
            isMixed = isMixed && scalarNear(result[i], ScalarType(2), 1E-10);
        expect(isMixed);
    }

    void testXCProvider() {
        const Index3D dim{4, 4, 4};
        LDAProvider lda(dim[0] * dim[1] * dim[2]);
        LDAProvider::DensityType density(1, dim);
        LDAProvider::PotType xc(1, dim);
        density[0].flatten() = ScalarType(1);

        lda.fill(density, xc);
        constexpr double exchangeFactor = -0.73855876638202240588;
        constexpr double correlationFactor1 = -0.045 / 2;
        constexpr double correlationFactor2 = 33.851831034345862;
        const ScalarType answer = ScalarType(4.0 / 3 * exchangeFactor)
                                + ScalarType(correlationFactor1) * ln(ScalarType(1) + ScalarType(correlationFactor2));
        const auto& potential = xc[0].flatten();
        bool isFilled = true;
        for (size_t i = 0; i < potential.getLength(); ++i)
            isFilled = isFilled && scalarNear(potential[i], answer, Precision);
        expect(isFilled);
    }

    void testKSSolver(const CrystalCellType& cell) {
        constexpr ScalarType cutEnergyPsi = 30;
        constexpr ScalarType cutEnergyRho = 200;
        const size_t numElectron = cell.getElectronCount();
        const size_t numBand = 2;
        BandType band(cell.makeRepLattice(), 1, 1, 1, numElectron, numBand);
        SolverType solver(cell, cutEnergyPsi, cutEnergyRho, band, numBand);

        expect(solver.getNumElectron() == numElectron);
        expect(solver.getNumFilledBand() == (numElectron + 1) / 2);
        expect(solver.getNumBand() == numBand);
        expect(solver.getNumPlainWave() > 0);

        SpinPair<BasisType, false> initial(cutEnergyPsi, cell.makeRepLattice());
        auto& coeff = initial[SpinState::Up].getCoeffGrid();
        forND(coeff.getShape(), [&](Index3D index) {
            coeff[index] = ComplexType(ScalarType(index[0]), ScalarType(2 * index[1]));
        });
        solver.initWaveFunc<RandomSource>(initial);

        const auto& orbitCoeff = solver.getOrbits()[0][SpinState::Up].getCoeffGrid();
        expect(orbitCoeff.getShape() == coeff.getShape());
        bool isCopied = true;
        forND(coeff.getShape(), [&](Index3D index) {
            isCopied = isCopied && orbitCoeff.calc(index) == coeff.calc(index);
        });
        expect(isCopied);
    }
}

int main() {
    const CrystalCellType cell = makeUnitCell();
    const LatticeMatrix repLatt = cell.makeRepLattice();

    testForKInBasis(repLatt);
    testDensityInterpolation();
    testBandGridSwap(repLatt);
    testKSHamilton(cell);
    testChargeMixer(repLatt);
    testXCProvider();
    testKSSolver(cell);
    return 0;
}
