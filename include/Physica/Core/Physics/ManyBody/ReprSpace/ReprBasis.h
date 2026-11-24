/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <type_traits>
#include "Physica/CRTP.h"

namespace Physica {
    template<class Derived>
    class ReprBasis : public CRTP<ReprBasis<Derived>> {
        using This = ReprBasis<Derived>;
        using Base = CRTP<This>;
    public:
        using StateType = Traits<Derived>::StateType;
        constexpr static unsigned int Dim = Traits<Derived>::Dim;
    public:
        ~ReprBasis() = default;
        /* Operators */
        [[nodiscard]] StateType operator[](size_t index) const noexcept { return Base::getDerived()[index]; }
        [[nodiscard]] size_t operator[](StateType state) const noexcept { return Base::getDerived()[state]; }
        /* Getters */
        [[nodiscard]] size_t getNumState() const noexcept { return Base::getDerived().getSize(); }
    protected:
        ReprBasis() = default;
        ReprBasis(const ReprBasis&) = default;
        ReprBasis(ReprBasis&&) noexcept = default;
        /* Operators */
        ReprBasis& operator=(const ReprBasis&) = default;
        ReprBasis& operator=(ReprBasis&&) noexcept = default;
    };

    template<class T>
    concept Representation = std::is_base_of<ReprBasis<T>, T>::value;
}

namespace Physica {
    template<class T>
    class Traits<ReprBasis<T>> {
    public:
        using Derived = T;
    };
}
