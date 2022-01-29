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

#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "BandGrid.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType, class XCProvider, bool isSpinPolarized> class KSSolverImpl;
    }

    template<class ScalarType, class XCProvider>
    class KSSolver : private Internal::KSSolverImpl<ScalarType, XCProvider, XCProvider::isSpinPolarized> {
        constexpr static bool isSpinPolarized = XCProvider::isSpinPolarized;
        using Base = Internal::KSSolverImpl<ScalarType, XCProvider, isSpinPolarized>;
        using typename Base::ComplexType;
        using typename Base::CenteredGrid;

        CenteredGrid externalPot;
    public:
        KSSolver(CrystalCell cell, ScalarType cutEnergy, BandGrid<ScalarType, isSpinPolarized> band, size_t gridDimX, size_t gridDimY, size_t gridDimZ);
        KSSolver(const KSSolver&) = delete;
        KSSolver(KSSolver&&) noexcept = delete;
        ~KSSolver() = default;
        /* Operators */
        KSSolver& operator=(const KSSolver& base) = delete;
        KSSolver& operator=(KSSolver&& base) noexcept = delete;
        /* Operations */
        bool solve(const ScalarType& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] const BandGrid<ScalarType, isSpinPolarized>& getBand() const noexcept { return Base::band; }
    private:
        /* Operations */
        void initExternalPot();
        /* Getters */
        [[nodiscard]] Utils::Array<CenteredGrid> getStructureFactor(ScalarType factorCutoff);
    };

    template<class ScalarType, class XCProvider>
    KSSolver<ScalarType, XCProvider>::KSSolver(CrystalCell cell,
                                               ScalarType cutEnergy,
                                               BandGrid<ScalarType, isSpinPolarized> band,
                                               size_t gridDimX,
                                               size_t gridDimY,
                                               size_t gridDimZ)
            : Base(std::move(cell), std::move(cutEnergy), std::move(band), gridDimX, gridDimY, gridDimZ) {
        initExternalPot();
    }

    template<class ScalarType, class XCProvider>
    bool KSSolver<ScalarType, XCProvider>::solve(const ScalarType& criteria, size_t maxIte) {
        return Base::solve(criteria, maxIte, externalPot);
    }

    template<class ScalarType, class XCProvider>
    void KSSolver<ScalarType, XCProvider>::initExternalPot() {
        const ScalarType factorCutoff = Base::cutEnergy * 8;
        externalPot = CenteredGrid::gridFromCutEnergy(factorCutoff, Base::repCell.getLattice());
        const Utils::Array<CenteredGrid> all_factors = getStructureFactor(factorCutoff);
        const ScalarType factor1 = ScalarType(-4 * M_PI) / Base::cell.getVolume();
        const std::unordered_set<uint16_t> species = Base::cell.getSpecies();

        externalPot.asVector() = ScalarType::Zero();
        const size_t gridSize = externalPot.getSize();

        size_t j = 0;
        for (uint16_t element : species) {
            const CenteredGrid& factors = all_factors[j];
            for (size_t i = 0; i < gridSize; ++i)
                externalPot[i] += factor1 * Base::getCharge(element) * factors[i] / factors.indexToPos(i).squaredNorm();
            ++j;
        }
        externalPot(0, 0, 0) = ComplexType::Zero();
    }

    template<class ScalarType, class XCProvider>
    Utils::Array<typename KSSolver<ScalarType, XCProvider>::CenteredGrid>
    KSSolver<ScalarType, XCProvider>::getStructureFactor(ScalarType factorCutoff) {
        const std::unordered_set<uint16_t> species = Base::cell.getSpecies();
        const auto& lattice = Base::repCell.getLattice();
        auto all_factors = Utils::Array<CenteredGrid>(species.size(), CenteredGrid::gridFromCutEnergy(factorCutoff, lattice));
        const size_t factors_size = all_factors[0].getSize();
        const size_t atomCount = Base::cell.getAtomCount();

        Vector<ScalarType, 3> g;
        size_t j = 0;
        for (uint16_t element : species) {
            CenteredGrid& factors = all_factors[j];
            for (size_t i = 0; i < factors_size; ++i) {
                g = factors.indexToPos(i);
                auto temp = ComplexType::Zero();
                for (size_t ion = 0; ion < atomCount; ++ion) {
                    if (Base::cell.getAtomicNumber(ion) == element) { //We can use searching table method
                        auto r = Base::cell.getPos().row(ion);
                        const ScalarType phase = g * r;
                        temp += ComplexType(cos(phase), sin(phase));
                    }
                }
                factors.asVector()[i] = temp;
            }
            ++j;
        }
        return all_factors;
    }
}

#include "KSSolverImpl.h"
