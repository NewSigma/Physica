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
     * [2] VASP (www.vasp.at)
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
        constexpr static double RSpaceSumPrec = 4;
        constexpr static double KSpaceSumPrec = 8;
    private:
        LatticeMatrix lattice;
        ReciprocalCell<PosScalarType> repCell;
        Vector<OutScalarType> charges;
        PosScalarType inv_volume;
        PosScalarType averageCellSize;
        OutScalarType integralLimit;
        OutScalarType selfEnergy;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
        Vector<OutScalarType> erfc_table;
    public:
        Ewald() = default;
        Ewald(LatticeMatrix lattice_, Vector<OutScalarType> charges_);
        Ewald(LatticeMatrix lattice_, ReciprocalCell<PosScalarType> repCell_, Vector<OutScalarType> charges_);
        Ewald(const Ewald&) = default;
        Ewald(Ewald&&) noexcept = default;
        ~Ewald() = default;
        /* Operators */
        Ewald& operator=(Ewald ewald) noexcept;
        [[nodiscard]] OutScalarType operator()(const PositionMatrix& pos) const;
        /* Operations */
        void swap(Ewald& ewald) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
    private:
        void makeErfcTable();
        OutScalarType erfcFromTable(OutScalarType x) const;
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
        /* The following param chosen is referenced from VASP, it is NOT the only choice */ {
            averageCellSize = cbrt(volume);
            integralLimit = sqrt(OutScalarType(M_PI)) / averageCellSize;
            rSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(lattice, OutScalarType(RSpaceSumPrec) / integralLimit);
            kSpaceSumRange = PeriodicCell<PosScalarType, Dim>::estimateRange(repCell.getLattice(), OutScalarType(KSpaceSumPrec) * integralLimit);
        }
        selfEnergy = square(charges).sum() * integralLimit / sqrt(OutScalarType(M_PI))
                   + square(charges.sum()) * OutScalarType(M_PI) / (OutScalarType::Two() * square(integralLimit)) * inv_volume;
        makeErfcTable();
    }

    template<class OutScalarType, class PosScalarType>
    Ewald<OutScalarType, PosScalarType>& Ewald<OutScalarType, PosScalarType>::operator=(Ewald ewald) noexcept {
        swap(ewald);
        return *this;
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::operator()(const PositionMatrix& pos) const {
        using ScalarType = OutScalarType;
        ScalarType rSpaceSum = 0;
        PeriodicCell<PosScalarType, Dim>::forParticleInRange(rSpaceSumRange, lattice, [this, pos, &rSpaceSum](Vector3D delta) {
                ScalarType sum = 0;
                for (size_t i = 0; i < getNumParticle(); ++i) {
                    const Vector3D pos_i = pos.row(i).asVector() + delta;
                    for (size_t j = i; j < getNumParticle(); ++j) {
                        const Vector3D pos_ij = pos_i - pos.row(j).asVector();
                        const ScalarType r2 = pos_ij.squaredNorm();
                        const bool isNotSelf = std::numeric_limits<ScalarType>::min() < r2;
                        if (isNotSelf) {
                            const ScalarType r = sqrt(r2);
                            const ScalarType temp = erfcFromTable(integralLimit * r) / r * (charges[i] * charges[j]); //Optimize: VASP uses searching table method
                            sum += temp * ScalarType(i == j ? 1 : 2);
                        }
                    }
                }
                rSpaceSum += sum;
            });
        
        ScalarType kSpaceSum = 0;
        PeriodicCell<PosScalarType, Dim>::forParticleInRange(kSpaceSumRange, repCell.getLattice(), [this, pos, &kSpaceSum](Vector3D delta) {
                const ScalarType squaredNorm = delta.squaredNorm();
                const bool isNotGammaPoint = std::numeric_limits<ScalarType>::min() < squaredNorm;
                if (isNotGammaPoint) {
                    ScalarType sum = 0;
                    for (size_t i = 0; i < getNumParticle(); ++i) {
                        for (size_t j = i; j < getNumParticle(); ++j) {
                            const ScalarType dot = delta * (pos.row(i).asVector() - pos.row(j).asVector());
                            const ScalarType temp = cos(dot) * (charges[i] * charges[j]);
                            sum += temp * ScalarType(i == j ? 1 : 2);
                        }
                    }
                    const ScalarType factor = reciprocal(squaredNorm * exp(squaredNorm / square(ScalarType::Two() * integralLimit)));
                    kSpaceSum += sum * factor;
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
        inv_volume.swap(ewald.inv_volume);
        averageCellSize.swap(ewald.averageCellSize);
        integralLimit.swap(ewald.integralLimit);
        selfEnergy.swap(ewald.selfEnergy);
        rSpaceSumRange.swap(ewald.rSpaceSumRange);
        kSpaceSumRange.swap(ewald.kSpaceSumRange);
        erfc_table.swap(ewald.erfc_table);
    }

    template<class OutScalarType, class PosScalarType>
    void Ewald<OutScalarType, PosScalarType>::makeErfcTable() {
        for (size_t i = 0; i < erfc_table.getLength(); ++i)
            erfc_table[i] = erfc(OutScalarType(i * TableStep));
    }

    template<class OutScalarType, class PosScalarType>
    OutScalarType Ewald<OutScalarType, PosScalarType>::erfcFromTable(OutScalarType x) const {
        using MatrixType = DenseMatrix<OutScalarType, MatrixOption::Row | MatrixOption::Element, Dim, Dim>;
        const size_t index = double(x * (1 / TableStep) + 0.5);
        if (index == 0)
            return 1;
        if (index + 1 >= TableSize)
            return 0;
        const OutScalarType x1 = TableStep * (index - 1);
        const OutScalarType x2 = TableStep * index;
        const OutScalarType x3 = TableStep * (index + 1);
        
        const MatrixType mat{square(x1), x1, 1, square(x2), x2, 1, square(x3), x3, 1};
        const Vector<OutScalarType, 3> y{erfc_table[index - 1], erfc_table[index], erfc_table[index + 1]};
        const Vector<OutScalarType, 3> coeff = MatrixType(mat.inverse()) * y;
        return (coeff[0] * x + coeff[1]) * x + coeff[2];
    }
}
