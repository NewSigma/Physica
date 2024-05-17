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
#include "State/SpinElectron.h"
#include "Physica/Utils/Container/Array/Array.h"
#include "Physica/Utils/CUDA/PlainStruct.h"

namespace Physica::Core {
    template<bool UseInversionSymm>
    class SpinRepr;

    namespace Internal {
        template<bool UseInversionSymm>
        class Traits<SpinRepr<UseInversionSymm>> {
        public:
            using StateType = SpinElectron;
            constexpr static unsigned int Dim = 1;
            constexpr static bool IsTransInvariant = false;
        };
    }

    template<bool UseInversionSymm>
    class SpinRepr : public ReprBasis<SpinRepr<UseInversionSymm>> {
        using This = SpinRepr<UseInversionSymm>;
        using Base = ReprBasis<This>;
    public:
        using typename Base::StateType;
    private:
        using StateArray = Utils::Array<SpinlessElectron>;
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
        [[nodiscard]] inline size_t getNumState() const noexcept;
    private:
        [[nodiscard]] StateArray makeSpinlessStates(size_t numElectron) const noexcept;
        void checkState(StateType state) const noexcept;
    };

    template<bool UseInversionSymm>
    SpinRepr<UseInversionSymm>::SpinRepr(unsigned int numSite_, unsigned int numSpinUp_, unsigned int numSpinDown_)
            : numSite(numSite_), numSpinUp(numSpinUp_), numSpinDown(numSpinDown_) {
        assert(!((numSpinUp == numSpinDown) ^ UseInversionSymm) && "[Error]: Inconsistent inversion symmetry");
        upStates = makeSpinlessStates(numSpinUp);
        if constexpr (!UseInversionSymm)
            downStates = makeSpinlessStates(numSpinDown);
    }

    template<bool UseInversionSymm>
    typename SpinRepr<UseInversionSymm>::StateType
    SpinRepr<UseInversionSymm>::operator[](size_t index) const noexcept {
        assert(index < getNumState() && "[Error]: Index out of range");
        const size_t upIndex = index / upStates.getLength();
        const size_t downIndex = index % upStates.getLength();
        StateType result;
        if constexpr (UseInversionSymm)
            result = StateType(upStates[upIndex], upStates[downIndex]);
        else
            result = StateType(upStates[upIndex], downStates[downIndex]);
        checkState(result);
        return result;
    }

    template<bool UseInversionSymm>
    size_t SpinRepr<UseInversionSymm>::operator[](StateType state) const noexcept {
        checkState(state);
        size_t upIndex = 0;
        for (; upIndex < upStates.getLength(); ++upIndex)
            if (state.getSpinUp() == upStates[upIndex])
                break;

        size_t downIndex = 0;
        if constexpr (UseInversionSymm) {
            for (; downIndex < upStates.getLength(); ++downIndex)
                if (state.getSpinDown() == upStates[downIndex])
                    break;
        }
        else {
            for (; downIndex < downStates.getLength(); ++downIndex)
                if (state.getSpinDown() == downStates[downIndex])
                    break;
        }

        const size_t index = upIndex * upStates.getLength() + downIndex;
        assert(index < getNumState() && "[Error]: Index out of range");
        return index;
    }

    template<bool UseInversionSymm>
    void SpinRepr<UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(numSpinUp, obj.numSpinUp);
        std::swap(numSpinDown, obj.numSpinDown);
        upStates.swap(obj.upStates);
        downStates.swap(obj.downStates);
    }

    template<bool UseInversionSymm>
    inline const typename SpinRepr<UseInversionSymm>::StateArray& SpinRepr<UseInversionSymm>::getDownStates() const noexcept {
        if constexpr (UseInversionSymm)
            return upStates;
        else
            return downStates;
    }

    template<bool UseInversionSymm>
    inline size_t SpinRepr<UseInversionSymm>::getNumState() const noexcept {
        const size_t numUpStates = upStates.getLength();
        if constexpr (UseInversionSymm)
            return numUpStates * numUpStates;
        else
            return numUpStates * downStates.getLength();
    }

    template<bool UseInversionSymm>
    typename SpinRepr<UseInversionSymm>::StateArray
    SpinRepr<UseInversionSymm>::makeSpinlessStates(size_t numElectron) const noexcept {
        const size_t numSpinlessState = SpinlessElectron::calcFullNumState(numSite);
        StateArray result{};
        result.reserve(numSpinlessState);
        for (size_t i = 0; i < numSpinlessState; ++i) {
            const SpinlessElectron state(i, numSite);
            if (state.getNumElectron() != numElectron)
                continue;
            result.append(state);
        }
        result.squeeze();
        return result;
    }

    template<bool UseInversionSymm>
    void SpinRepr<UseInversionSymm>::checkState([[maybe_unused]] StateType state) const noexcept {
        assert(state.getNumSpinUpElectron() == numSpinUp && "[Error]: Unexpected state");
        assert(state.getNumSpinDownElectron() == numSpinDown && "[Error]: Unexpected state");
    }
}
