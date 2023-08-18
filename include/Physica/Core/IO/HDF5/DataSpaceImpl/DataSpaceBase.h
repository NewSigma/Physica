/*
 * Copyright 2023 WeiBo He.
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

    template<class Derived>
    class DataSpaceBase : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        constexpr static size_t Dim = Internal::Traits<Derived>::Dim;
        using SizeArray = Utils::Array<hsize_t, Dim>;
    public:
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) { Base::getDerived().selectHyperslab(op, count, start); }
        /* Getters */
        [[nodiscard]] const H5::DataSpace& asH5Type() const noexcept { return Base::getDerived().asH5Type(); }
        [[nodiscard]] H5::DataSpace& asH5Type() noexcept { return Base::getDerived().asH5Type(); }
        [[nodiscard]] size_t getDim() const noexcept { return Base::getDerived().getDim(); }
        [[nodiscard]] size_t getSize(size_t dim) const noexcept { return Base::getDerived().getSize(dim); }
        [[nodiscard]] const SizeArray& getSelectedCount() const noexcept { return Base::getDerived().getSelectedCount(); }
        [[nodiscard]] const SizeArray& getSelectedStart() const noexcept { return Base::getDerived().getSelectedStart(); }
    };
}
