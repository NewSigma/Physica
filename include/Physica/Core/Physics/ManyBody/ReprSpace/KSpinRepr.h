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
    template<bool UseInversionSymm> class KSpinRepr;

    namespace Internal {
        template<bool UseInversionSymm>
        class Traits<KSpinRepr<UseInversionSymm>> : public SpinRepr<UseInversionSymm> {
        public:
            constexpr static bool IsTransInvariant = true;
        };
    }

    template<bool UseInversionSymm>
    class KSpinRepr : public ReprBasis<KSpinRepr<UseInversionSymm>> {
        using This = KSpinRepr<UseInversionSymm>;
        using Base = ReprBasis<This>;
        using RSpinType = SpinRepr<UseInversionSymm>;
        using PeriodArray = Utils::Array<int>;
    public:
        using typename Base::StateType;
    private:
        Utils::Array<SpinElectron> states;
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
        [[nodiscard]] size_t operator[](StateType state) const noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] unsigned int getNumSite() const noexcept { return rSpin.getNumSite(); }
        [[nodiscard]] size_t getNumState() const noexcept { return states.getLength(); }
        [[nodiscard]] const PeriodArray& getPeriods() const noexcept { return periods; }
        [[nodiscard]] unsigned int getKIndex() const noexcept { return kIndex; }
        [[nodiscard]] inline unsigned int getReducedK() const noexcept;
    };

    template<bool UseInversionSymm>
    KSpinRepr<UseInversionSymm>::KSpinRepr(RSpinType rSpin_, unsigned int kIndex_) : rSpin(std::move(rSpin_)), kIndex(kIndex_) {
        assert(kIndex < getNumSite() && "[Error]: Momentum index out of range");
        const size_t numSite = getNumSite();
        for (auto psiUp : rSpin.getUpStates()) {
            if (psiUp.isTransReducible())
                continue;

            const int periodUp = psiUp.calcPeriod();
            for (auto psiDown : rSpin.getDownStates()) {
                if (psiDown.isTransReducible(periodUp))
                    continue;
                const int period = periodUp * psiDown.calcPeriod();
                const bool isZero = (kIndex * period) % numSite != 0;
                if (isZero)
                    continue;
                states.append(SpinElectron(psiUp, psiDown));
                periods.append(period);
            }
        }
    }

    template<bool UseInversionSymm>
    size_t KSpinRepr<UseInversionSymm>::operator[](StateType state) const noexcept {
        for (size_t i = 0; i < getNumState(); ++i)
            if (state == states[i])
                return i;
        assert(false && "[Error]: Invalid state");
        return 0;
    }

    template<bool UseInversionSymm>
    void KSpinRepr<UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        states.swap(obj.states);
        periods.swap(obj.periods);
        rSpin.swap(obj.rSpin);
        std::swap(kIndex, obj.kIndex);
    }

    template<bool UseInversionSymm>
    inline unsigned int KSpinRepr<UseInversionSymm>::getReducedK() const noexcept {
        const auto kSize = FFT<Scalar<>, 1>::rSizeToKSize(getNumSite());
        if (kIndex < kSize)
            return kIndex;
        return kSize - kIndex;
    }
}
