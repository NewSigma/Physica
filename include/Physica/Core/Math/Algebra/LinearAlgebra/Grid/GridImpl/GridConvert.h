/*
 * Copyright 2023-2024 Weibo He.
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

namespace Physica::Core {
    template<class GridType> class RValueGrid;

    template<class T>
    class RealGrid : public RValueGrid<RealGrid<T>> {
        using Base = RValueGrid<RealGrid<T>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit RealGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).real(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T>
    class ImagGrid : public RValueGrid<ImagGrid<T>> {
        using Base = RValueGrid<ImagGrid<T>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit ImagGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).imag(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T>
    class NormGrid : public RValueGrid<NormGrid<T>> {
        using Base = RValueGrid<NormGrid<T>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit NormGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).norm(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T>
    class ValueGrid : public RValueGrid<ValueGrid<T>> {
        using Base = RValueGrid<ValueGrid<T>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit ValueGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).value(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class T, int GradOrder>
    class GradGrid : public RValueGrid<GradGrid<T, GradOrder>> {
        using Base = RValueGrid<GradGrid<T, GradOrder>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit GradGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).template grad<GradOrder>(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };
}

namespace Physica {
    template<class T>
    class Traits<RealGrid<T>> {
    public:
        using ScalarType = T::ScalarType::RealType;
    };

    template<class T>
    class Traits<ImagGrid<T>> : public Traits<RealGrid<T>> {};

    template<class T>
    class Traits<NormGrid<T>> : public Traits<RealGrid<T>> {};

    template<class T>
    class Traits<ValueGrid<T>> {
        static_assert(T::ScalarType::isDiffable, "[Error]: Unnecessary toValueGrid() call or toGradGrid() call");
    public:
        using ScalarType = T::PlainType;
    };

    template<class T, int GradOrder>
    class Traits<GradGrid<T, GradOrder>> : public Traits<ValueGrid<T>> {};
}
