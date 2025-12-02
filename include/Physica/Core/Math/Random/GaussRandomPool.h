/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Scalar T, class RandomPoolType>
    class GaussRandomPool {
    public:
        using Generator = RandomPoolType::Generator;
    private:
        VectorND<T> rands;
    public:
        GaussRandomPool(size_t size);
        GaussRandomPool(const GaussRandomPool&) = default;
        GaussRandomPool(GaussRandomPool&&) noexcept = default;
        ~GaussRandomPool() = default;
        /* Operators */
        GaussRandomPool& operator=(GaussRandomPool obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T operator()() const noexcept;
        /* Operations */
        template<RNG R>
        void init();
        void swap(GaussRandomPool& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return rands.getLength(); }
        /* Static member */
        [[nodiscard]] static Generator& generator() noexcept { return RandomPoolType::getInstance().generator(); }
    };

    template<Scalar T, class RandomPoolType>
    GaussRandomPool<T, RandomPoolType>::GaussRandomPool(size_t size) : rands(size) {
        init(RandomPoolType::getInstance().generator());
    }

    template<Scalar T, class RandomPoolType>
    template<RNG R>
    void GaussRandomPool<T, RandomPoolType>::init() {
        rands.template random_normal<R>();
    }

    template<Scalar T, class RandomPoolType>
    T GaussRandomPool<T, RandomPoolType>::operator()() const noexcept {
        std::uniform_int_distribution<size_t> dist(0, getSize() - 1);
        return rands[dist(RandomPoolType::getInstance().generator())];
    }

    template<Scalar T, class RandomPoolType>
    void GaussRandomPool<T, RandomPoolType>::swap(GaussRandomPool& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        rands.swap(obj.rands);
    }
}

namespace Physica {
    template<Scalar T, class RandomPoolType>
    class Traits<GaussRandomPool<T, RandomPoolType>> : public Traits<RandomPoolType> {};
}
