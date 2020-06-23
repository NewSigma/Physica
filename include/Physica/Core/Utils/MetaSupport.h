/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_METASUPPORT_H
#define PHYSICA_METASUPPORT_H

namespace Physica::Core {
    template<size_t s1, size_t s2>
    struct META_MAX
    { static constexpr size_t value = s1 > s2 ? s1 : s2; };

    template<size_t s1, size_t s2>
    struct META_MIN
    { static constexpr size_t value = s1 > s2 ? s2 : s1; };
}

#endif