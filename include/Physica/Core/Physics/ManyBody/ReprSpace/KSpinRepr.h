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

#include "SpinRepr.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica::Core {
    template<unsigned int Dim, bool UseInversionSymm> class KSpinRepr;

    namespace Internal {
        template<unsigned int Dim, bool UseInversionSymm>
        class Traits<KSpinRepr<Dim, UseInversionSymm>> : public SpinRepr<Dim, UseInversionSymm> {
        public:
            constexpr static bool IsTransInvariant = true;
        };
    }

    template<unsigned int Dim, bool UseInversionSymm>
    class KSpinRepr : public ReprBasis<KSpinRepr<Dim, UseInversionSymm>> {
        using This = KSpinRepr<Dim, UseInversionSymm>;
        using Base = ReprBasis<This>;
        using RSpinType = SpinRepr<Dim, UseInversionSymm>;
        using PeriodArray = Utils::Array<int>;
    public:
        using typename Base::StateType;
    private:
        Utils::Array<StateType> states; //Optimize: Each spin up state might pair with several spin down states
        PeriodArray periods;
        RSpinType rSpin;
        unsigned int kIndex;
    public:
        KSpinRepr() = default;
        KSpinRepr(RSpinType rSpin_, unsigned int kIndex_);
        KSpinRepr(const KSpinRepr&) = default;
        KSpinRepr(KSpinRepr&&) noexcept = default;
        ~KSpinRepr() = default;
        /* Operators */
        KSpinRepr& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] StateType operator[](size_t index) const noexcept { return states[index]; }
        [[nodiscard]] size_t operator[](StateType psi) const noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] unsigned int getNumSite() const noexcept { return rSpin.getNumSite(); }
        [[nodiscard]] size_t getNumState() const noexcept { return states.getLength(); }
        [[nodiscard]] const PeriodArray& getPeriods() const noexcept { return periods; }
        [[nodiscard]] unsigned int getKIndex() const noexcept { return kIndex; }
        [[nodiscard]] inline unsigned int getReducedK() const noexcept;
    };

    template<unsigned int Dim, bool UseInversionSymm>
    KSpinRepr<Dim, UseInversionSymm>::KSpinRepr(RSpinType rSpin_, unsigned int kIndex_) : rSpin(std::move(rSpin_)), kIndex(kIndex_) {
        assert(kIndex < getNumSite() && "[Error]: Momentum index out of range");
        const size_t numSite = getNumSite();
        for (auto psiUp : rSpin.getUpStates()) {
            if (psiUp.isTransReducible())
                continue;

            const int periodUp = psiUp.calcPeriod();
            for (auto psiDown : rSpin.getDownStates()) {
                if (psiDown.isTransReducible(periodUp))
                    continue;
                const int period = lcm<int, false>(periodUp, psiDown.calcPeriod());
                const bool isZero = (kIndex * period) % numSite != 0;
                if (isZero)
                    continue;
                states.append(StateType(psiUp, psiDown));
                periods.append(period);
            }
        }
    }

    template<unsigned int Dim, bool UseInversionSymm>
    size_t KSpinRepr<Dim, UseInversionSymm>::operator[](StateType psi) const noexcept {
        size_t left = 0, right = getNumState() - 1;  
        while (left < right) {
            const size_t mid = left + (right - left) / 2;
            const auto& psi0 = states[mid];
            if (psi0 == psi)
                return mid;
            if (psi0 < psi)
                left = mid + 1;
            else
                right = mid - 1;
        }
        return left;
    }

    template<unsigned int Dim, bool UseInversionSymm>
    void KSpinRepr<Dim, UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        states.swap(obj.states);
        periods.swap(obj.periods);
        rSpin.swap(obj.rSpin);
        std::swap(kIndex, obj.kIndex);
    }

    template<unsigned int Dim, bool UseInversionSymm>
    inline unsigned int KSpinRepr<Dim, UseInversionSymm>::getReducedK() const noexcept {
        const auto kSize = FFT<Scalar<>, 1>::rSizeToKSize(getNumSite());
        if (kIndex < kSize)
            return kIndex;
        return getNumSite() - kIndex;
    }
}
