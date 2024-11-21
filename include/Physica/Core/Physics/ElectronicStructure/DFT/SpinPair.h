/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Utils/Container/Array.h"

namespace Physica::Core {
    enum class SpinState {
        Up = 0,
        Down = 1
    };

    template<class T, bool IsSpinPolarized>
    class SpinPair : public Traits<SpinPair<T, IsSpinPolarized>>::Base {
        using Base = Traits<SpinPair<T, IsSpinPolarized>>::Base;
        constexpr static unsigned char NumSpin = Base::getLength();
    public:
        using ElemType = T;
    public:
        template<class... Args>
        explicit SpinPair(Args... args); //TODO: Use Args&&... may lead to incorrect overload
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
    template<class... Args>
    SpinPair<T, IsSpinPolarized>::SpinPair(Args... args) : Base(NumSpin, std::forward<Args>(args)...) {}
}

namespace Physica {
    template<class T, bool IsSpinPolarized>
    class Traits<Core::SpinPair<T, IsSpinPolarized>> {
    public:
        constexpr static unsigned char NumSpin = IsSpinPolarized ? 2 : 1;
        using Base = Array<T, NumSpin>;
    };
}
