/*
 * Copyright 2024-2026 Weibo He.
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

#include <concepts>
#include <type_traits>

namespace Physica {
    namespace Internal {
        template<template<class...> class, class>
        struct is_instance_of : std::false_type {};

        template<template<class...> class Template, class... Args>
        struct is_instance_of<Template, Template<Args...>> : std::true_type {};

        template<template<auto, class...> class, class>
        struct is_instance_of_x : std::false_type {};

        template<template<auto, class...> class Template, auto Arg0, class... Args>
        struct is_instance_of_x<Template, Template<Arg0, Args...>> : std::true_type {};

        template<template<auto, auto, class...> class, class>
        struct is_instance_of_xx : std::false_type {};

        template<template<auto, auto, class...> class Template, auto Arg0, auto Arg1, class... Args>
        struct is_instance_of_xx<Template, Template<Arg0, Arg1, Args...>> : std::true_type {};

        template<template<class, auto...> class, class>
        struct is_instance_of_tx : std::false_type {};

        template<template<class, auto...> class Template, class Arg0, auto... Args>
        struct is_instance_of_tx<Template, Template<Arg0, Args...>> : std::true_type {};

        template<template<class, class, auto...> class, class>
        struct is_instance_of_ttx : std::false_type {};

        template<template<class, class, auto...> class Template, class Arg0, class Arg1, auto... Args>
        struct is_instance_of_ttx<Template, Template<Arg0, Arg1, Args...>> : std::true_type {};

        template<template<class, auto, auto, auto, class...> class, class>
        struct is_instance_of_txxx : std::false_type {};

        template<template<class, auto, auto, auto, class...> class Template, class Arg0, auto Arg1, auto Arg2, auto Arg3, class... Args>
        struct is_instance_of_txxx<Template, Template<Arg0, Arg1, Arg2, Arg3, Args...>> : std::true_type {};
    }

    template<class T, template<class...> class Template>
    concept instanceof = Internal::is_instance_of<Template, std::remove_cvref_t<T>>::value;
    /**
     * x: Non type template param
     * t: Type template param
     *
     * It is unfortunate that we cannot implement instanceof in an elegant way.
     */
    template<class T, template<auto, class...> class Template>
    concept instanceof_x = Internal::is_instance_of_x<Template, std::remove_cvref_t<T>>::value;

    template<class T, template<auto, auto, class...> class Template>
    concept instanceof_xx = Internal::is_instance_of_xx<Template, std::remove_cvref_t<T>>::value;

    template<class T, template<class, auto...> class Template>
    concept instanceof_tx = Internal::is_instance_of_tx<Template, std::remove_cvref_t<T>>::value;

    template<class T, template<class, class, auto...> class Template>
    concept instanceof_ttx = Internal::is_instance_of_ttx<Template, std::remove_cvref_t<T>>::value;

    template<class T, template<class, auto, auto, auto, class...> class Template>
    concept instanceof_txxx = Internal::is_instance_of_txxx<Template, std::remove_cvref_t<T>>::value;
    /**
     * Reject const&& to avoid potential bad pattern
     */
    template<class Ref> requires(std::is_lvalue_reference_v<Ref> || (std::is_rvalue_reference_v<Ref> && !std::is_const_v<std::remove_reference_t<Ref>>))
    using LazyDestroy = std::conditional<std::is_rvalue_reference<Ref>::value, std::remove_reference_t<Ref>, Ref>::type;
    /**
     * std::forward_like does not work if we use GCC 14.2.
     * FIXME: Simply remove it in the future
     *
     * Reference:
     * [1] GH101614; https://bgithub.xyz/llvm/llvm-project/issues/101614
     */
    template<typename Tp, typename Up>
    [[nodiscard, gnu::nodebug]] constexpr decltype(auto) forward_like(Up&& x) noexcept {
        using namespace std;
        constexpr bool _as_rval = is_rvalue_reference_v<Tp&&>;
        if constexpr (is_const_v<remove_reference_t<Tp>>) {
            using Up2 = remove_reference_t<Up>;
            if constexpr (_as_rval)
                return static_cast<const Up2&&>(x);
            else
                return static_cast<const Up2&>(x);
        }
        else {
            if constexpr (_as_rval)
                return static_cast<remove_reference_t<Up>&&>(x);
            else
                return x;
        }
    }
    /**
     * Consider:
     * 1. A view that owns data
     * 2. We would like to move the data and create a new view
     * This function returns correct type if the move is eligible.
     *
     * \returns a rvalue reference if both is rvalue reference, otherwise a lvalue reference
     */
    template<class T1, class T2>
    [[nodiscard, gnu::nodebug]] constexpr decltype(auto) propagate_rvalue_reference(auto&& x) noexcept {
        static_assert(std::same_as<std::remove_cvref_t<T2>, std::remove_cvref_t<decltype(x)>>);
        static_assert(std::is_reference<T1>::value);
        static_assert(std::is_reference<T2>::value);
        if constexpr (std::is_rvalue_reference_v<T1> && std::is_rvalue_reference_v<T2>)
            return forward_like<T2>(x);
        else
            return x;
    }

    template<class From, class To>
    class copy_cvref {
        using T0 = std::remove_cvref<To>::type;
        using T1 = std::conditional<std::is_const_v<std::remove_reference_t<From>>, const T0, T0>::type;
        using T2 = std::conditional<std::is_rvalue_reference_v<From>, T1&&, T1>::type;
    public:
        using type = std::conditional<std::is_lvalue_reference_v<From>, T2&, T2>::type;
    };
}
