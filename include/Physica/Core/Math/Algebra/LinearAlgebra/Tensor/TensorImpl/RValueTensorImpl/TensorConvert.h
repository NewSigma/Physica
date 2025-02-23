/*
 * Copyright 2023-2025 Weibo He.
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

namespace Physica {
    template<class TensorType> class RValueTensor;

    template<class T>
    class RealTensor : public RValueTensor<RealTensor<T>> {
        using Base = RValueTensor<RealTensor<T>>;
        const T& g;
    public:
        using typename Base::ScalarType;
    public:
        explicit RealTensor(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).real(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T>
    class ImagTensor : public RValueTensor<ImagTensor<T>> {
        using Base = RValueTensor<ImagTensor<T>>;
        const T& g;
    public:
        using typename Base::ScalarType;
    public:
        explicit ImagTensor(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).imag(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T>
    class NormTensor : public RValueTensor<NormTensor<T>> {
        using Base = RValueTensor<NormTensor<T>>;
        const T& g;
    public:
        using typename Base::ScalarType;
    public:
        explicit NormTensor(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).norm(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T>
    class ValueTensor : public RValueTensor<ValueTensor<T>> {
        using Base = RValueTensor<ValueTensor<T>>;
        const T& g;
    public:
        using typename Base::ScalarType;
    public:
        explicit ValueTensor(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).value(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T, int GradOrder>
    class GradTensor : public RValueTensor<GradTensor<T, GradOrder>> {
        using Base = RValueTensor<GradTensor<T, GradOrder>>;
        const T& g;
    public:
        using typename Base::ScalarType;
    public:
        explicit GradTensor(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).template grad<GradOrder>(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<RealTensor<T>> {
    public:
        using ScalarType = T::ScalarType::RealType;
    };

    template<class T>
    class Traits<ImagTensor<T>> : public Traits<RealTensor<T>> {};

    template<class T>
    class Traits<NormTensor<T>> : public Traits<RealTensor<T>> {};

    template<class T>
    class Traits<ValueTensor<T>> {
        static_assert(T::ScalarType::isDiffable, "[Error]: Unnecessary toValueTensor() call or toGradTensor() call");
    public:
        using ScalarType = T::PlainType;
    };

    template<class T, int GradOrder>
    class Traits<GradTensor<T, GradOrder>> : public Traits<ValueTensor<T>> {};
}
