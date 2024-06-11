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

#include "ReprImpl/ReprSpace.h"
#include "State/SpinFermion.h"
#include "Physica/Utils/Container/Array/Array.h"
#include "Physica/Utils/CUDA/PlainStruct.h"

namespace Physica::Core {
    template<unsigned int Dim, bool UseInversionSymm>
    class SpinRepr;

    namespace Internal {
        template<unsigned int I, bool UseInversionSymm>
        class Traits<SpinRepr<I, UseInversionSymm>> {
        public:
            using StateType = SpinFermion<I>;
            constexpr static unsigned int Dim = I;
            constexpr static bool IsTransInvariant = false;
        };
    }

    template<unsigned int Dim, bool UseInversionSymm>
    class SpinRepr : public ReprBasis<SpinRepr<Dim, UseInversionSymm>> {
        using This = SpinRepr<Dim, UseInversionSymm>;
        using Base = ReprBasis<This>;
    public:
        using typename Base::StateType;
    private:
        using SpinlessState = SpinlessFermion<Dim>;
        using StateArray = Utils::Array<SpinlessState>;
        using DownStateArray = typename std::conditional<UseInversionSymm, PlainStruct<void>, StateArray>::type;

        unsigned int numSite;
        unsigned int numSpinUp;
        unsigned int numSpinDown;
        StateArray upStates;
        DownStateArray downStates;
    public:
        SpinRepr() = default;
        SpinRepr(unsigned int numSite_, unsigned int numSpinUp_, unsigned int numSpinDown_);
        SpinRepr(const SpinRepr&) = default;
        SpinRepr(SpinRepr&&) noexcept = default;
        ~SpinRepr() = default;
        /* Operators */
        SpinRepr& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] StateType operator[](size_t index) const noexcept;
        [[nodiscard]] size_t operator[](StateType state) const noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] unsigned int getNumSite() const noexcept { return numSite; }
        [[nodiscard]] const StateArray& getUpStates() const noexcept { return upStates; }
        [[nodiscard]] inline const StateArray& getDownStates() const noexcept;
        [[nodiscard]] size_t getNumState() const noexcept { return getNumUpStates() * getNumDownStates(); }
        [[nodiscard]] size_t getNumUpStates() const noexcept { return getUpStates().getLength(); }
        [[nodiscard]] size_t getNumDownStates() const noexcept { return getDownStates().getLength(); }
    private:
        [[nodiscard]] StateArray makeSpinlessStates(size_t numElectron) const noexcept;
        void checkState(StateType state) const noexcept;
    };

    template<unsigned int Dim, bool UseInversionSymm>
    SpinRepr<Dim, UseInversionSymm>::SpinRepr(unsigned int numSite_, unsigned int numSpinUp_, unsigned int numSpinDown_)
            : numSite(numSite_), numSpinUp(numSpinUp_), numSpinDown(numSpinDown_) {
        assert(!((numSpinUp == numSpinDown) ^ UseInversionSymm) && "[Error]: Inconsistent inversion symmetry");
        upStates = makeSpinlessStates(numSpinUp);
        if constexpr (!UseInversionSymm)
            downStates = makeSpinlessStates(numSpinDown);
    }

    template<unsigned int Dim, bool UseInversionSymm>
    typename SpinRepr<Dim, UseInversionSymm>::StateType
    SpinRepr<Dim, UseInversionSymm>::operator[](size_t index) const noexcept {
        assert(index < getNumState() && "[Error]: Index out of range");
        const size_t numDownStates = getNumDownStates();
        const size_t upIndex = index / numDownStates;
        const size_t downIndex = index % numDownStates;
        StateType result = StateType(upStates[upIndex], getDownStates()[downIndex]);
        checkState(result);
        return result;
    }

    template<unsigned int Dim, bool UseInversionSymm>
    size_t SpinRepr<Dim, UseInversionSymm>::operator[](StateType state) const noexcept {
        checkState(state);
        size_t upIndex = 0;
        for (; upIndex < upStates.getLength(); ++upIndex)
            if (state.getSpinUp() == upStates[upIndex])
                break;
        assert(upIndex < upStates.getLength() && "[Error]: Unexpected missing state");

        const size_t numDownStates = getNumDownStates();
        size_t downIndex = 0;
        for (; downIndex < numDownStates; ++downIndex)
            if (state.getSpinDown() == getDownStates()[downIndex])
                break;
        assert(downIndex < numDownStates && "[Error]: Unexpected missing state");

        const size_t index = upIndex * getNumDownStates() + downIndex;
        assert(index < getNumState() && "[Error]: Index out of range");
        return index;
    }

    template<unsigned int Dim, bool UseInversionSymm>
    void SpinRepr<Dim, UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(numSpinUp, obj.numSpinUp);
        std::swap(numSpinDown, obj.numSpinDown);
        upStates.swap(obj.upStates);
        downStates.swap(obj.downStates);
    }

    template<unsigned int Dim, bool UseInversionSymm>
    inline const typename SpinRepr<Dim, UseInversionSymm>::StateArray& SpinRepr<Dim, UseInversionSymm>::getDownStates() const noexcept {
        if constexpr (UseInversionSymm)
            return upStates;
        else
            return downStates;
    }

    template<unsigned int Dim, bool UseInversionSymm>
    typename SpinRepr<Dim, UseInversionSymm>::StateArray
    SpinRepr<Dim, UseInversionSymm>::makeSpinlessStates(size_t numElectron) const noexcept {
        const size_t numSpinlessState = SpinlessState::calcFullNumState(numSite);
        StateArray result{};
        result.reserve(numSpinlessState);
        for (size_t i = 0; i < numSpinlessState; ++i) {
            const SpinlessState state(i, numSite);
            if (state.getNumElectron() != numElectron)
                continue;
            result.append(state);
        }
        result.squeeze();
        return result;
    }

    template<unsigned int Dim, bool UseInversionSymm>
    void SpinRepr<Dim, UseInversionSymm>::checkState([[maybe_unused]] StateType state) const noexcept {
        assert(state.getNumSpinUpElectron() == numSpinUp && "[Error]: Unexpected state");
        assert(state.getNumSpinDownElectron() == numSpinDown && "[Error]: Unexpected state");
    }
}
