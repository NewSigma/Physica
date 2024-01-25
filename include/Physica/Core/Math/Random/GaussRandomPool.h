/*
 * Copyright 2024 WeiBo He.
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

#include <random>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class ScalarType, class RandomPoolType> class GaussRandomPool;

    namespace Internal {
        template<class ScalarType, class RandomPoolType>
        class Traits<GaussRandomPool<ScalarType, RandomPoolType>> : public Traits<RandomPoolType> {};
    }

    template<class ScalarType, class RandomPoolType>
    class GaussRandomPool {
        Vector<ScalarType> rands;
    public:
        GaussRandomPool(size_t size);
        GaussRandomPool(const GaussRandomPool&) = default;
        GaussRandomPool(GaussRandomPool&&) noexcept = default;
        ~GaussRandomPool() = default;
        /* Operators */
        GaussRandomPool& operator=(GaussRandomPool obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] inline ScalarType operator()() const noexcept;
        /* Operations */
        template<class RandomGenerator>
        inline void init(RandomGenerator& gen);
        inline void swap(GaussRandomPool& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return rands.getLength(); }
    };

    template<class ScalarType, class RandomPoolType>
    GaussRandomPool<ScalarType, RandomPoolType>::GaussRandomPool(size_t size) : rands(size) {
        init(RandomPoolType::getInstance().getGen());
    }

    template<class ScalarType, class RandomPoolType>
    template<class RandomGenerator>
    inline void GaussRandomPool<ScalarType, RandomPoolType>::init(RandomGenerator& gen) {
        rands.random_normal(gen);
    }

    template<class ScalarType, class RandomPoolType>
    inline ScalarType GaussRandomPool<ScalarType, RandomPoolType>::operator()() const noexcept {
        std::uniform_int_distribution<size_t> dist(0, getSize() - 1);
        return rands[dist(RandomPoolType::getInstance().getGen())];
    }

    template<class ScalarType, class RandomPoolType>
    inline void GaussRandomPool<ScalarType, RandomPoolType>::swap(GaussRandomPool& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        rands.swap(obj.rands);
    }
}
