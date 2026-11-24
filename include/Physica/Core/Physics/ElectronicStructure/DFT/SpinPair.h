/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    enum class SpinState : int8_t {
        Up = 0,
        Down = 1
    };

    template<class T, bool IsSpinPolarized>
    class SpinPair : public Traits<SpinPair<T, IsSpinPolarized>>::Base {
        using Base = Traits<SpinPair<T, IsSpinPolarized>>::Base;
        constexpr static unsigned char NumSpin = Base::getLength();
    public:
        using value_type = T;
    public:
        explicit SpinPair(auto... args); //TODO: Use auto&&... may lead to incorrect overload
        SpinPair(const SpinPair&) = default;
        SpinPair(SpinPair&&) noexcept = default;
        ~SpinPair() = default;
        /* Operators */
        SpinPair& operator=(SpinPair pair) noexcept { swap(pair); return *this; }
        [[nodiscard]] T& operator[](SpinState spin) { return Base::operator[](IsSpinPolarized ? int(spin) : 0); }
        [[nodiscard]] const T& operator[](SpinState spin) const { return Base::operator[](IsSpinPolarized ? int(spin) : 0); }
        /* Operations */
        void swap(SpinPair& __restrict pair) noexcept { Base::swap(pair); }
    private:
        using Base::operator[];
    };

    template<class T, bool IsSpinPolarized>
    SpinPair<T, IsSpinPolarized>::SpinPair(auto... args) : Base(NumSpin, std::forward<decltype(args)>(args)...) {}
}

namespace Physica {
    template<class T, bool IsSpinPolarized>
    class Traits<SpinPair<T, IsSpinPolarized>> {
    public:
        constexpr static unsigned char NumSpin = IsSpinPolarized ? 2 : 1;
        using Base = Array<T, NumSpin>;
    };
}
