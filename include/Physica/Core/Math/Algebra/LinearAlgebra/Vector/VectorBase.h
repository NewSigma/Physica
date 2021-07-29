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
    /**
     * \class VectorBase is the parent class of a collection of classes that can be assigned to a vector.
     */
    template<class Derived>
    class VectorBase : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
    public:
        template<class OtherDerived>
        void assignTo(VectorBase<OtherDerived>& v) const {
            assert(v.getLength() == getLength());
            Base::getDerived().assignTo(v);
        }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
    };
}
