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
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/PlainStruct.h"
#include "State/FermiState.h"
#include "ReprBasis.h"

namespace Physica {
    template<int Dim, int NumSite, bool UseInversionSymm>
    class FermiRepr : public ReprBasis<FermiRepr<Dim, NumSite, UseInversionSymm>> {
        using This = FermiRepr<Dim, NumSite, UseInversionSymm>;
        using Base = ReprBasis<This>;
    public:
        using typename Base::StateType;
    private:
        using Spin = SpinState<Dim, NumSite>;
        using StateArray = Array<Spin>;
        using DownStateArray = std::conditional<UseInversionSymm, PlainStruct<void>, StateArray>::type;
        using StateIndexMap = std::unordered_map<Spin, size_t>;
        using DownStateIndexMap = std::conditional<UseInversionSymm, PlainStruct<void>, StateIndexMap>::type;
        using PairType = std::pair<StateArray, StateIndexMap>;

        int numSpinUp;
        int numSpinDown;
        StateArray upStates;
        StateIndexMap upIndexMap;
        [[no_unique_address]] DownStateArray downStates;
        [[no_unique_address]] DownStateIndexMap downIndexMap;
    public:
        FermiRepr() = default;
        FermiRepr(int numSpinUp_, int numSpinDown_);
        FermiRepr(const This&) = default;
        FermiRepr(This&&) noexcept = default;
        ~FermiRepr() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] StateType operator[](size_t index) const noexcept;
        [[nodiscard]] size_t operator[](StateType state) const noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getNumSpinUp() const noexcept { return numSpinUp; }
        [[nodiscard]] int getNumSpinDown() const noexcept { return numSpinDown; }
        [[nodiscard]] int getNumParticle() const noexcept { return numSpinUp + numSpinDown; }
        [[nodiscard]] const auto& getUpStates() const noexcept { return upStates; }
        [[nodiscard]] const auto& getUpIndexMap() const noexcept { return upIndexMap; }
        [[nodiscard]] const StateArray& getDownStates() const noexcept;
        [[nodiscard]] const StateIndexMap& getDownIndexMap() const noexcept;
        [[nodiscard]] size_t getNumState() const noexcept { return getNumUpStates() * getNumDownStates(); }
        [[nodiscard]] size_t getNumUpStates() const noexcept { return getUpStates().getLength(); }
        [[nodiscard]] size_t getNumDownStates() const noexcept { return getDownStates().getLength(); }
    private:
        [[nodiscard]] PairType makeSpinStates(int numParticle) const noexcept;
        void checkState(StateType state) const noexcept;
    };

    template<int Dim, int NumSite, bool UseInversionSymm>
    FermiRepr<Dim, NumSite, UseInversionSymm>::FermiRepr(int numSpinUp_, int numSpinDown_)
            : numSpinUp(numSpinUp_), numSpinDown(numSpinDown_) {
        assert((numSpinUp >= 0) && (numSpinDown >= 0));
        assert(((numSpinUp == numSpinDown) == UseInversionSymm) && "[Error]: Inconsistent inversion symmetry");
        auto pair = makeSpinStates(numSpinUp);
        upStates = std::move(pair.first);
        upIndexMap = std::move(pair.second);
        if constexpr (!UseInversionSymm) {
            auto pair = makeSpinStates(numSpinDown);
            downStates = std::move(pair.first);
            downIndexMap = std::move(pair.second);
        }
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    auto FermiRepr<Dim, NumSite, UseInversionSymm>::operator[](size_t index) const noexcept -> StateType {
        assert(index < getNumState() && "[Error]: Index out of range");
        const size_t numDownStates = getNumDownStates();
        const size_t upIndex = index / numDownStates;
        const size_t downIndex = index % numDownStates;
        StateType result = StateType(upStates[upIndex], getDownStates()[downIndex]);
        checkState(result);
        return result;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    size_t FermiRepr<Dim, NumSite, UseInversionSymm>::operator[](StateType state) const noexcept {
        checkState(state);
        const size_t upIndex = upIndexMap.find(state.getSpinUp())->second;
        const size_t downIndex = getDownIndexMap().find(state.getSpinDown())->second;
        const size_t index = upIndex * getNumDownStates() + downIndex;
        assert(index < getNumState() && "[Error]: Index out of range");
        return index;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    void FermiRepr<Dim, NumSite, UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(numSpinUp, obj.numSpinUp);
        std::swap(numSpinDown, obj.numSpinDown);
        upStates.swap(obj.upStates);
        downStates.swap(obj.downStates);
        upIndexMap.swap(obj.upIndexMap);
        downIndexMap.swap(obj.downIndexMap);
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    auto FermiRepr<Dim, NumSite, UseInversionSymm>::getDownStates() const noexcept -> const StateArray& {
        if constexpr (UseInversionSymm)
            return upStates;
        else
            return downStates;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    auto FermiRepr<Dim, NumSite, UseInversionSymm>::getDownIndexMap() const noexcept -> const StateIndexMap& {
        if constexpr (UseInversionSymm)
            return getUpIndexMap();
        else
            return downIndexMap;
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    auto FermiRepr<Dim, NumSite, UseInversionSymm>::makeSpinStates(int numParticle) const noexcept -> PairType {
        constexpr size_t numSpinState = Spin::calcFullNumState();
        assert(numParticle >= 0);
        StateArray arr{};
        StateIndexMap map{};
        arr.reserve(numSpinState);
        for (size_t i = 0; i < numSpinState; ++i) {
            const Spin state(i);
            if (state.getNumParticle() != numParticle)
                continue;
            map[state] = arr.getLength();
            arr.append(state);
        }
        arr.squeeze();
        return std::make_pair(std::move(arr), std::move(map));
    }

    template<int Dim, int NumSite, bool UseInversionSymm>
    void FermiRepr<Dim, NumSite, UseInversionSymm>::checkState([[maybe_unused]] StateType state) const noexcept {
        assert(state.getNumSpinUpParticle() == numSpinUp && "[Error]: Unexpected state");
        assert(state.getNumSpinDownParticle() == numSpinDown && "[Error]: Unexpected state");
    }
}

namespace Physica {
    template<int I1, int I2, bool UseInversionSymm>
    class Traits<FermiRepr<I1, I2, UseInversionSymm>> {
    public:
        using StateType = FermiState<I1, I2>;
        constexpr static int Dim = I1;
        constexpr static bool IsTransInvariant = false;
    };
}
