/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Transform/FFT.h"
#include "FermiRepr.h"

namespace Physica {
    template<unsigned int Dim, unsigned int NumSite, bool UseInversionSymm>
    class KFermiRepr : public ReprBasis<KFermiRepr<Dim, NumSite, UseInversionSymm>> {
        using This = KFermiRepr<Dim, NumSite, UseInversionSymm>;
        using Base = ReprBasis<This>;
        using RSpinType = FermiRepr<Dim, NumSite, UseInversionSymm>;
        using PeriodArray = Array<int>;
    public:
        using typename Base::StateType;
    private:
        Array<StateType> states; //Optimize: Each spin up state might pair with several spin down states
        PeriodArray periods;
        RSpinType rSpin;
        unsigned int kIndex;
    public:
        KFermiRepr() = default;
        KFermiRepr(RSpinType rSpin_, unsigned int kIndex_);
        KFermiRepr(const This&) = default;
        KFermiRepr(This&&) noexcept = default;
        ~KFermiRepr() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] StateType operator[](size_t index) const noexcept { return states[index]; }
        [[nodiscard]] size_t operator[](StateType psi) const noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumState() const noexcept { return states.getLength(); }
        [[nodiscard]] const PeriodArray& getPeriods() const noexcept { return periods; }
        [[nodiscard]] unsigned int getKIndex() const noexcept { return kIndex; }
        [[nodiscard]] unsigned int getReducedK() const noexcept;
    };

    template<unsigned int Dim, unsigned int NumSite, bool UseInversionSymm>
    KFermiRepr<Dim, NumSite, UseInversionSymm>::KFermiRepr(RSpinType rSpin_, unsigned int kIndex_)
            : rSpin(std::move(rSpin_)), kIndex(kIndex_) {
        assert(kIndex < NumSite && "[Error]: Momentum index out of range");
        const size_t numSite = NumSite;
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

    template<unsigned int Dim, unsigned int NumSite, bool UseInversionSymm>
    size_t KFermiRepr<Dim, NumSite, UseInversionSymm>::operator[](StateType psi) const noexcept {
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

    template<unsigned int Dim, unsigned int NumSite, bool UseInversionSymm>
    void KFermiRepr<Dim, NumSite, UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        states.swap(obj.states);
        periods.swap(obj.periods);
        rSpin.swap(obj.rSpin);
        std::swap(kIndex, obj.kIndex);
    }

    template<unsigned int Dim, unsigned int NumSite, bool UseInversionSymm>
    unsigned int KFermiRepr<Dim, NumSite, UseInversionSymm>::getReducedK() const noexcept {
        const auto kSize = FFT<Real<>, 1>::rSizeToKSize(NumSite);
        if (kIndex < kSize)
            return kIndex;
        return NumSite - kIndex;
    }
}

namespace Physica {
    template<unsigned int Dim, unsigned int NumSite, bool UseInversionSymm>
    class Traits<KFermiRepr<Dim, NumSite, UseInversionSymm>> : public FermiRepr<Dim, NumSite, UseInversionSymm> {
    public:
        constexpr static bool IsTransInvariant = true;
    };
}
