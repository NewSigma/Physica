/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Grid3D.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:499-503
     * [2] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications; San Diego: Academic, 2002:304-306
     * [3] VASP (www.vasp.at)
     */
    template<class OutScalarType, class PosScalarType = OutScalarType>
    class Ewald {
        constexpr static unsigned int Dim = 3;
        using LatticeMatrix = typename PeriodicCell<PosScalarType, Dim>::LatticeMatrix;
        using PositionMatrix = typename PeriodicCell<PosScalarType, Dim>::PositionMatrix;
        using SearchRangeType = typename PeriodicCell<PosScalarType, Dim>::SearchRangeType;
        using Vector3D = Vector<PosScalarType, Dim>;
        constexpr static size_t TableSize = 4608;
        constexpr static double TableStep = 0.001;
        constexpr static double SumPrec = 4; // Refer to [2]
    private:
        LatticeMatrix lattice;
        ReciprocalCell<PosScalarType> repCell;
        Vector<OutScalarType> charges;
        OutScalarType sumSquaredCharges;
        PosScalarType inv_volume;
        OutScalarType integralLimit;
        OutScalarType selfEnergy;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
        Vector<OutScalarType> erfc_table;
        OutScalarType step;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice_, Vector<OutScalarType> charges_);
        Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<OutScalarType> charges_);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald ewald) noexcept;
        [[nodiscard]] Vector<OutScalarType> force(const PositionMatrix& pos) const;
        [[nodiscard]] OutScalarType operator()(const PositionMatrix& pos) const;
        /* Operations */
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
    private:
        void makeModifiedErfcTable();
        OutScalarType modifiedErfcFromTable(OutScalarType x) const;
    };

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, Vector<OutScalarType> charges_)
            : Ewald(lattice_, ReciprocalCell(lattice_), std::move(charges_)) {}

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>::Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<OutScalarType> charges_)
            : lattice(std::move(lattice_))
            , repCell(std::move(repCell_))
            , charges(std::move(charges_))
            , erfc_table(TableSize) {
        const PosScalarType volume = PeriodicCell<PosScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(volume);
        {
            const OutScalarType averageCellSize = cbrt(volume);
            integralLimit = sqrt(OutScalarType(M_PI)) / averageCellSize;
            rSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(lattice, PosScalarType(OutScalarType(SumPrec) / integralLimit));
            kSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(repCell.getLattice(), PosScalarType(OutScalarType(SumPrec * 2) * integralLimit));
        }
        sumSquaredCharges = square(charges).sum();
        selfEnergy = sumSquaredCharges * integralLimit / sqrt(OutScalarType(M_PI))
                   + square(charges.sum()) * OutScalarType(M_PI) / (OutScalarType::Two() * square(integralLimit)) * inv_volume;
        makeModifiedErfcTable();
    }

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>& Ewald<OutScalarType, PosScalarType>::operator=(Ewald ewald) noexcept {
        swap(ewald);
        return *this;
    }

    template<class OutScalarType, class PosScalarType>
    Vector<OutScalarType> Ewald<OutScalarType, PosScalarType>::force(const PositionMatrix& pos) const {
        using ScalarType = OutScalarType;
        const size_t numParticle = getNumParticle();
        Vector<OutScalarType> result(numParticle * Dim, 0);
        
        const OutScalarType factor1 = reciprocal(square(ScalarType::Two() * integralLimit));
        PeriodicCell<PosScalarType, Dim>::forParticleInRange(kSpaceSumRange, repCell.getLattice(), [this, numParticle, factor1, &pos, &result](Vector3D delta) {
                const ScalarType squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = std::numeric_limits<ScalarType>::min() < squaredNorm;
                if (isNotGammaPoint) {
                    const ScalarType factor2 = reciprocal(squaredNorm * exp(squaredNorm * factor1));
                    for (size_t i = 0; i < numParticle - 1; ++i) {
                        auto force_i = result.segment(i * Dim, (i + 1) * Dim);
                        const ScalarType charge_i = charges[i];
                        const Vector3D pos_i = pos.row(i).asVector();
                        for (size_t j = i + 1; j < numParticle; ++j) {
                            auto force_j = result.segment(j * Dim, (j + 1) * Dim);
                            const ScalarType dot = delta * (pos_i - pos.row(j).asVector());
                            const Vector<ScalarType, Dim> f = (sin(dot) * charge_i * charges[j] * factor2) * delta;
                            force_i += f;
                            force_j -= f;
                        }
                    }
                }
            });
        result *= ScalarType(4 * M_PI) * inv_volume;

        PeriodicCell<PosScalarType, Dim>::forParticleInRange(rSpaceSumRange, lattice, [this, numParticle, &pos, &result](Vector3D delta) {
                const ScalarType factor3 = ScalarType(2) / sqrt(ScalarType(M_PI)) * integralLimit;
                for (size_t i = 0; i < numParticle; ++i) {
                    const Vector3D pos_i = pos.row(i).asVector() + delta;
                    Vector<ScalarType, Dim> sum(Dim, 0);
                    for (size_t j = 0; j < numParticle; ++j) {
                        if (j == i)
                            continue;
                        const Vector3D pos_ij = pos_i - pos.row(j).asVector();
                        const ScalarType r2 = pos_ij.squaredNorm();
                        const ScalarType r = sqrt(r2);
                        const ScalarType temp = modifiedErfcFromTable(r) + factor3 * exp(-square(integralLimit * r));
                        sum += (temp / r2 * charges[j]) * pos_ij;
                    }
                    auto force_i = result.segment(i * Dim, (i + 1) * Dim);
                    force_i += charges[i] * sum;
                }
            });
        return result;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::operator()(const PositionMatrix& pos) const {
        using ScalarType = OutScalarType;
        const size_t numParticle = getNumParticle();

        ScalarType rSpaceSum = 0;
        PeriodicCell<PosScalarType, Dim>::forParticleInRange(rSpaceSumRange, lattice, [this, numParticle, &pos, &rSpaceSum](Vector3D delta) {
                ScalarType sum = 0;
                for (size_t i = 0; i < numParticle; ++i) {
                    const Vector3D pos_i = pos.row(i).asVector() + delta;
                    const OutScalarType charge_i = charges[i];
                    for (size_t j = i; j < numParticle; ++j) {
                        const Vector3D pos_ij = pos_i - pos.row(j).asVector();
                        const ScalarType r2 = pos_ij.squaredNorm();
                        const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                        if (isNotSelf) {
                            const ScalarType r = sqrt(r2);
                            const ScalarType temp = modifiedErfcFromTable(r) * (charge_i * charges[j]);
                            sum += temp * ScalarType(i == j ? 1 : 2);
                        }
                    }
                }
                rSpaceSum += sum;
            });
        
        ScalarType kSpaceSum = 0;
        const OutScalarType factor = reciprocal(square(ScalarType::Two() * integralLimit));
        PeriodicCell<PosScalarType, Dim>::forParticleInRange(kSpaceSumRange, repCell.getLattice(), [this, numParticle, factor, &pos, &kSpaceSum](Vector3D delta) {
                const ScalarType squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = std::numeric_limits<ScalarType>::min() < squaredNorm;
                if (isNotGammaPoint) {
                    ScalarType sum = 0;
                    for (size_t i = 0; i < numParticle - 1; ++i) {
                        const Vector3D pos_i = pos.row(i).asVector();
                        const OutScalarType charge_i = charges[i];
                        for (size_t j = i + 1; j < numParticle; ++j) {
                            const ScalarType dot = delta * (pos_i - pos.row(j).asVector());
                            const ScalarType temp = cos(dot) * (charge_i * charges[j]);
                            sum += temp;
                        }
                    }
                    kSpaceSum += (sumSquaredCharges + sum * 2) / (squaredNorm * exp(squaredNorm * factor));
                }
            });
        kSpaceSum *= ScalarType(4 * M_PI) * inv_volume;
        return (rSpaceSum + kSpaceSum) * 0.5 - selfEnergy;
    }

    template<class OutScalarType, class PosScalarType>
    void Ewald<OutScalarType, PosScalarType>::swap(Ewald& ewald) noexcept {
        lattice.swap(ewald.lattice);
        repCell.swap(ewald.repCell);
        charges.swap(ewald.charges);
        sumSquaredCharges.swap(ewald.sumSquaredCharges);
        inv_volume.swap(ewald.inv_volume);
        integralLimit.swap(ewald.integralLimit);
        selfEnergy.swap(ewald.selfEnergy);
        rSpaceSumRange.swap(ewald.rSpaceSumRange);
        kSpaceSumRange.swap(ewald.kSpaceSumRange);
        erfc_table.swap(ewald.erfc_table);
    }

    template<class OutScalarType, class PosScalarType>
    void Ewald<OutScalarType, PosScalarType>::makeModifiedErfcTable() {
        for (size_t i = 0; i < erfc_table.getLength(); ++i) {
            OutScalarType x = OutScalarType(i * TableStep);
            erfc_table[i] = erfc(x) / x * integralLimit;
        }
        step = OutScalarType(TableStep) / integralLimit;
    }

    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::modifiedErfcFromTable(OutScalarType x) const {
        using MatrixType = DenseMatrix<OutScalarType, MatrixOption::Row | MatrixOption::Element, Dim, Dim>;
        const size_t index = double(x / step + 0.5);
        if (index == 0)
            return 1;
        if (index + 1 >= TableSize)
            return 0;
        const OutScalarType x2 = step * index;
        const OutScalarType x1 = x2 - step;
        const OutScalarType x3 = x2 + step;
        
        const MatrixType mat{square(x1), x1, 1, square(x2), x2, 1, square(x3), x3, 1};
        const Vector<OutScalarType, 3> y = erfc_table.segment(index - 1, index + 2);
        const Vector<OutScalarType, 3> coeff = MatrixType(mat.inverse()) * y;
        return ((coeff[0] * x + coeff[1]) * x + coeff[2]);
    }
}
