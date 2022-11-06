/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    template<class T, bool isSpinPolarized> class SpinPair;

    namespace Internal {
        template<class T, bool isSpinPolarized>
        class Traits<SpinPair<T, isSpinPolarized>> {
        public:
            constexpr static unsigned char NumSpin = isSpinPolarized ? 2 : 1;
            using Base = Utils::Array<T, NumSpin>;
        };
    }

    enum class SpinState {
        Up = 0,
        Down = 1
    };

    template<class T, bool isSpinPolarized>
    class SpinPair : public Internal::Traits<SpinPair<T, isSpinPolarized>>::Base {
        using Base = typename Internal::Traits<SpinPair<T, isSpinPolarized>>::Base;
        constexpr static unsigned char NumSpin = Base::getLength();
    public:
        using ElemType = T;
    public:
        template<class... Args>
        explicit SpinPair(Args... args);
        SpinPair(const SpinPair&) = default;
        SpinPair(SpinPair&&) noexcept = default;
        ~SpinPair() = default;
        /* Operators */
        SpinPair& operator=(SpinPair pair) noexcept { swap(pair); return *this; }
        [[nodiscard]] T& operator[](SpinState spin) { return Base::operator[](isSpinPolarized ? int(spin) : 0); }
        [[nodiscard]] const T& operator[](SpinState spin) const { return Base::operator[](isSpinPolarized ? int(spin) : 0); }
        /* Helpers */
        void swap(SpinPair& pair) noexcept { Base::swap(pair); }
    private:
        using Base::operator[];
    };

    template<class T, bool isSpinPolarized>
    template<class... Args>
    SpinPair<T, isSpinPolarized>::SpinPair(Args... args) : Base(NumSpin, std::forward<Args>(args)...) {}
}
