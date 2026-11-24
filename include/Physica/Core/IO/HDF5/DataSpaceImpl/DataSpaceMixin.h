/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <array>
#include <cassert>
#include "Physica/CRTP.h"
#include <H5Spublic.h>

namespace Physica {
    template<class Derived>
    class DataSpaceMixin : public CRTP<DataSpaceMixin<Derived>> {
        using This = DataSpaceMixin<Derived>;
        using Base = CRTP<This>;
    public:
        constexpr static size_t Dim = Traits<Derived>::Dim;
        using SizeType = hsize_t;
        using SizeArray = std::array<hsize_t, Dim>;
    public:
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) { Base::getDerived().selectHyperslab(op, count, start); }
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return Base::getDerived().getDim(); }
        [[nodiscard]] size_t getSize(size_t dim) const noexcept { return Base::getDerived().getSize(dim); }
        [[nodiscard]] const SizeArray& getSelectedCount() const noexcept { return Base::getDerived().getSelectedCount(); }
        [[nodiscard]] const SizeArray& getSelectedStart() const noexcept { return Base::getDerived().getSelectedStart(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<DataSpaceMixin<T>> {
    public:
        using Derived = T;
    };
}
