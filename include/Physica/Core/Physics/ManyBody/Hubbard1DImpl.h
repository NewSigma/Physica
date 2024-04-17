/*
 * Copyright 2024 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType>
    Hubbard1D<ScalarType>::Hubbard1D(LatticeType lattice, ScalarType hoppingT_, ScalarType repelU_, unsigned int numSpinUp_, unsigned int numSpinDown_)
            : Base(std::move(lattice)), hoppingT(hoppingT_), repelU(repelU_), numSpinUp(numSpinUp_), numSpinDown(numSpinDown_) {}

    template<class ScalarType>
    template<class VectorType>
    Vector<ScalarType> Hubbard1D<ScalarType>::operator*(const RValueVector<VectorType>& v) const {
        const size_t length = v.getLength();
        const auto numSite = getNumSuperCellSite();
        assert(Base::getColumn() == length && "[Error]: Dimensions do not match");
        Vector<ScalarType> result(length, 0);
        for (size_t i = 0; i < length; ++i) {
            const ScalarType hoppingElem = v.calc(i) * hoppingT;
            const auto state = indexToState(i);
            int numRepel = 0;
            for (unsigned int site = 0; site < numSite; ++site) {
                const auto site1 = (site + 1) % numSite;
                stateAdd(result, state.hopUp(site, site1), hoppingElem);
                stateAdd(result, state.hopUp(site1, site), hoppingElem);
                stateAdd(result, state.hopDown(site, site1), hoppingElem);
                stateAdd(result, state.hopDown(site1, site), hoppingElem);
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            }
            result[i] += v.calc(i) * repelU * ScalarType(numRepel);
        }
        return result;
    }

    template<class ScalarType>
    ScalarType Hubbard1D<ScalarType>::calc(size_t row, size_t col) const {
        if (row == 0 || col == 0) //TODO: Implement particle num converve and remove this
            return ScalarType(0);

        const auto colState = indexToState(col);
        const auto numSite = getNumSuperCellSite();
        if (row == col) {
            int numRepel = 0;
            for (unsigned int site = 0; site < numSite; ++site)
                numRepel += colState.isUpOccupy(site) && colState.isDownOccupy(site);
            return repelU * ScalarType(numRepel);
        }

        const auto rowState = indexToState(row);
        for (unsigned int site = 0; site < numSite; ++site) {
            const auto site1 = (site + 1) % numSite;
            if (rowState == colState.hopUp(site, site1)
             || rowState == colState.hopUp(site1, site)
             || rowState == colState.hopDown(site, site1)
             || rowState == colState.hopDown(site1, site))
                return -hoppingT;
        }
        return ScalarType(0);
    }

    template<class ScalarType>
    void Hubbard1D<ScalarType>::swap(Hubbard1D& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        std::swap(numSpinUp, obj.numSpinUp);
        std::swap(numSpinDown, obj.numSpinDown);
    }

    template<class ScalarType>
    size_t Hubbard1D<ScalarType>::stateToIndex(StateType state) const noexcept {
        checkState(state);
        const size_t index = state.getSpinUp().getOccupyBits() | (state.getSpinDown().getOccupyBits() << getNumSuperCellSite());
        assert(index < Base::getNumState() && "[Error]: Index out of range");
        return index;
    }

    template<class ScalarType>
    typename Hubbard1D<ScalarType>::StateType Hubbard1D<ScalarType>::indexToState(size_t index) const noexcept {
        assert(index < Base::getNumState() && "[Error]: Index out of range");
        const size_t spinDownBits = index >> getNumSuperCellSite();
        const size_t spinUpBits = index ^ (spinDownBits << getNumSuperCellSite());
        StateType result(spinUpBits, spinDownBits);
        checkState(result);
        return result;
    }

    template<class ScalarType>
    void Hubbard1D<ScalarType>::stateAdd(Vector<ScalarType>& target, StateType state, ScalarType value) const noexcept {
        if (state.isVacuum())
            return;
        const auto index = stateToIndex(state);
        target[index] += value;
    }

    template<class ScalarType>
    void Hubbard1D<ScalarType>::checkState([[maybe_unused]] StateType state) const noexcept {
        //assert(state.getNumSpinUpElectron() == numSpinUp && "[Error]: Unexpected state");
        //assert(state.getNumSpinDownElectron() == numSpinDown && "[Error]: Unexpected state");
    }
}
