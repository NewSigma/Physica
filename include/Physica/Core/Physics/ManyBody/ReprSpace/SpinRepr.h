/*
 * Copyright 2024 Weibo He.
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

#include <unordered_map>
#include "ReprImpl/ReprSpace.h"
#include "State/SpinFermion.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/PlainStruct.h"

namespace Physica::Core {
    template<int Dim, int NumSite, bool UseInversionSymm>
    class SpinRepr : public ReprBasis<SpinRepr<Dim, NumSite, UseInversionSymm>> {
        using This = SpinRepr<Dim, NumSite, UseInversionSymm>;
        using Base = ReprBasis<This>;
    public:
        using typename Base::StateType;
    private:
        using SpinlessState = SpinlessFermion<Dim, NumSite>;
        using StateArray = Array<SpinlessState>;
        using DownStateArray = typename std::conditional<UseInversionSymm, PlainStruct<void>, StateArray>::type;
        using StateIndexMap = std::unordered_map<SpinlessState, size_t>;
        using DownStateIndexMap = typename std::conditional<UseInversionSymm, PlainStruct<void>, StateIndexMap>::type;
        using PairType = std::pair<StateArray, StateIndexMap>;

        int numSpinUp;
        int numSpinDown;
        StateArray upStates;
        DownStateArray downStates;
        StateIndexMap upIndexMap;
        DownStateIndexMap downIndexMap;
    public:
        SpinRepr() = default;
        SpinRepr(int numSpinUp_, int numSpinDown_);
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
        [[nodiscard]] int getNumSpinUp() const noexcept { return numSpinUp; }
        [[nodiscard]] int getNumSpinDown() const noexcept { return numSpinDown; }
        [[nodiscard]] int getNumParticle() const noexcept { return numSpinUp + numSpinDown; }
        [[nodiscard]] const StateArray& getUpStates() const noexcept { return upStates; }
        [[nodiscard]] inline const StateArray& getDownStates() const noexcept;
        [[nodiscard]] const StateIndexMap& getUpIndexMap() const noexcept { return upIndexMap; }
        [[nodiscard]] inline const StateIndexMap& getDownIndexMap() const noexcept;
        [[nodiscard]] size_t getNumState() const noexcept { return getNumUpStates() * getNumDownStates(); }
        [[nodiscard]] size_t getNumUpStates() const noexcept { return getUpStates().getLength(); }
        [[nodiscard]] size_t getNumDownStates() const noexcept { return getDownStates().getLength(); }
    private:
        [[nodiscard]] PairType makeSpinlessStates(int numParticle) const noexcept;
        [[nodiscard]] size_t findStateIndex(const StateArray& arr, SpinlessState psi) const;
        void checkState(StateType state) const noexcept;
    };

    template<int Dim, int NumSite, bool UseInversionSymm>
    SpinRepr<Dim, NumSite, UseInversionSymm>::SpinRepr(int numSpinUp_, int numSpinDown_)
            : numSpinUp(numSpinUp_), numSpinDown(numSpinDown_) {
        assert((numSpinUp >= 0) && (numSpinDown >= 0));
        assert(!((numSpinUp == numSpinDown) ^ UseInversionSymm) && "[Error]: Inconsistent inversion symmetry");
        auto pair = makeSpinlessStates(numSpinUp);
        upStates = std::move(pair.first);
        upIndexMap = std::move(pair.second);
        if constexpr (!UseInversionSymm) {
            auto pair = makeSpinlessStates(numSpinDown);
            downStates = std::move(pair.first);
            downIndexMap = std::move(pair.second);
        }
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    typename SpinRepr<Dim, NumSite, UseInversionSymm>::StateType
    SpinRepr<Dim, NumSite, UseInversionSymm>::operator[](size_t index) const noexcept {
        assert(index < getNumState() && "[Error]: Index out of range");
        const size_t numDownStates = getNumDownStates();
        const size_t upIndex = index / numDownStates;
        const size_t downIndex = index % numDownStates;
        StateType result = StateType(upStates[upIndex], getDownStates()[downIndex]);
        checkState(result);
        return result;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    size_t SpinRepr<Dim, NumSite, UseInversionSymm>::operator[](StateType state) const noexcept {
        checkState(state);
        const size_t upIndex = upIndexMap.find(state.getSpinUp())->second;
        const size_t downIndex = getDownIndexMap().find(state.getSpinDown())->second;
        const size_t index = upIndex * getNumDownStates() + downIndex;
        assert(index < getNumState() && "[Error]: Index out of range");
        return index;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    void SpinRepr<Dim, NumSite, UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(numSpinUp, obj.numSpinUp);
        std::swap(numSpinDown, obj.numSpinDown);
        upStates.swap(obj.upStates);
        downStates.swap(obj.downStates);
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    inline const typename SpinRepr<Dim, NumSite, UseInversionSymm>::StateArray&
    SpinRepr<Dim, NumSite, UseInversionSymm>::getDownStates() const noexcept {
        if constexpr (UseInversionSymm)
            return upStates;
        else
            return downStates;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    inline const typename SpinRepr<Dim, NumSite, UseInversionSymm>::StateIndexMap&
    SpinRepr<Dim, NumSite, UseInversionSymm>::getDownIndexMap() const noexcept {
        if constexpr (UseInversionSymm)
            return getUpIndexMap();
        else
            return downIndexMap;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    typename SpinRepr<Dim, NumSite, UseInversionSymm>::PairType
    SpinRepr<Dim, NumSite, UseInversionSymm>::makeSpinlessStates(int numParticle) const noexcept {
        constexpr size_t numSpinlessState = SpinlessState::calcFullNumState();
        assert(numParticle >= 0);
        StateArray arr{};
        StateIndexMap map{};
        arr.reserve(numSpinlessState);
        for (size_t i = 0; i < numSpinlessState; ++i) {
            const SpinlessState state(i);
            if (state.getNumParticle() != numParticle)
                continue;
            map[state] = arr.getLength();
            arr.append(state);
        }
        arr.squeeze();
        return std::make_pair(std::move(arr), std::move(map));
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    size_t SpinRepr<Dim, NumSite, UseInversionSymm>::findStateIndex(const StateArray& arr, SpinlessState psi) const {
        size_t left = 0, right = arr.getLength() - 1;  
        while (left < right) {
            const size_t mid = left + (right - left) / 2;
            const auto& psi0 = arr[mid];
            if (psi0 == psi)
                return mid;
            if (psi0 < psi)
                left = mid + 1;
            else
                right = mid - 1;
        }
        assert(psi == arr[left] && "[Error]: Unexpected missing state");
        return left;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    void SpinRepr<Dim, NumSite, UseInversionSymm>::checkState([[maybe_unused]] StateType state) const noexcept {
        assert(state.getNumSpinUpParticle() == numSpinUp && "[Error]: Unexpected state");
        assert(state.getNumSpinDownParticle() == numSpinDown && "[Error]: Unexpected state");
    }
}

namespace Physica {
    template<int I1, int I2, bool UseInversionSymm>
    class Traits<Core::SpinRepr<I1, I2, UseInversionSymm>> {
    public:
        using StateType = Core::SpinFermion<I1, I2>;
        constexpr static int Dim = I1;
        constexpr static bool IsTransInvariant = false;
    };
}
