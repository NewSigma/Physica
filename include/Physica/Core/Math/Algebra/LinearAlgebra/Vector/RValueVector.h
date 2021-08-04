/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    namespace Internal {
        template<class T> class Traits;
    }

    template<class Derived> class LValueVector;
    /**
     * \class RValueVector is base class of vectors that can be assigned to \class LValueVector
     * but other vectors cannot be assigned to this class.
     * In other words, you cannot take the address of elements in the vector.
     */
    template<class Derived>
    class RValueVector : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
    public:
        template<class OtherDerived>
        void assignTo(LValueVector<OtherDerived>& v) const {
            assert(v.getLength() == getLength());
            Base::getDerived().assignTo(v);
        }

        [[nodiscard]] ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
    };
}
