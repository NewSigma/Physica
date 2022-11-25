/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Physics/PIPhonon.h"

namespace Physica::Core {
    PIPhonon::PIPhonon(size_t numAtomUnitCell_,
                       size_t superSizeX_,
                       size_t superSizeY_,
                       size_t superSizeZ_)
            : fft({superSizeX_, superSizeY_, superSizeZ_}, {1, 1, 1})
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
    }

    PIPhonon& PIPhonon::operator=(PIPhonon obj) noexcept {
        swap(obj);
        return *this;
    }

    std::ostream& operator<<(std::ostream& os, const PIPhonon& phonon) {
        for (size_t r = 0; r < phonon.getUnitCellDOF(); ++r) {
            for (size_t c = r; c < phonon.getUnitCellDOF(); ++c) {
                const size_t offset = phonon.force_corr.accessingIndex(r, c);
                os << phonon.force_corr[offset] << phonon.momentum_corr[offset];
            }
        }
        return os;
    }

    std::istream& operator>>(std::istream& is, PIPhonon& phonon) {
        for (size_t r = 0; r < phonon.getUnitCellDOF(); ++r) {
            for (size_t c = r; c < phonon.getUnitCellDOF(); ++c) {
                const size_t offset = phonon.force_corr.accessingIndex(r, c);
                is >> phonon.force_corr[offset] >> phonon.momentum_corr[offset];
            }
        }
        return is;
    }

    void PIPhonon::swap(PIPhonon& obj) noexcept {
        fft.swap(obj.fft);
        std::swap(numAtomUnitCell, obj.numAtomUnitCell);
        std::swap(superSizeX, obj.superSizeX);
        std::swap(superSizeY, obj.superSizeY);
        std::swap(superSizeZ, obj.superSizeZ);
        force_corr.swap(obj.force_corr);
        momentum_corr.swap(obj.momentum_corr);
        std::swap(numSample, obj.numSample);
    }

    void PIPhonon::toKSpace() {
        const size_t dof = getUnitCellDOF();
        for (size_t r = 0; r < dof; ++r) {
            for (size_t c = r; c < dof; ++c) {
                const size_t offset_corr = force_corr.accessingIndex(r, c);
                auto* corr = &force_corr[offset_corr];
                fft.transform(*corr);
                //*corr = toRealVector(fft.getFreqs());

                corr = &momentum_corr[offset_corr];
                fft.transform(*corr);
                //*corr = toRealVector(fft.getFreqs());
            }
        }
    }
}
