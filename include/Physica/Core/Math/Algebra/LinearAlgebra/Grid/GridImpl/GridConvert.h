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

    template<Grid T>
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

    template<Grid T>
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

    template<Grid T>
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

    template<Grid T>
    class ValueGrid : public RValueGrid<ValueGrid<T>> {
        using Base = RValueGrid<ValueGrid<T>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit ValueGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).getValue(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<Grid T>
    class GradGrid : public RValueGrid<GradGrid<T>> {
        using Base = RValueGrid<GradGrid<T>>;
        const T& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit GradGrid(const T& g_) : g(g_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).getGrad(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<Grid T>
    [[nodiscard]] inline RealGrid<T> toRealGrid(const T& v) {
        return RealGrid<T>{v};
    }

    template<Grid T>
    [[nodiscard]] inline ImagGrid<T> toImagGrid(const T& v) {
        return ImagGrid<T>{v};
    }

    template<Grid T>
    [[nodiscard]] inline NormGrid<T> toNormGrid(const T& v) {
        return NormGrid<T>{v};
    }

    template<Grid T>
    [[nodiscard]] inline ValueGrid<T> toValueGrid(const T& v) {
        return ValueGrid<T>{v};
    }

    template<Grid T>
    [[nodiscard]] inline GradGrid<T> toGradGrid(const T& v) {
        return GradGrid<T>{v};
    }
}

namespace Physica {
    template<Grid T>
    class Traits<RealGrid<T>> {
    public:
        using ScalarType = T::ScalarType::RealType;
    };

    template<Grid T>
    class Traits<ImagGrid<T>> : public Traits<RealGrid<T>> {};

    template<Grid T>
    class Traits<NormGrid<T>> : public Traits<RealGrid<T>> {};

    template<Grid T>
    class Traits<ValueGrid<T>> {
        static_assert(T::ScalarType::isDifferentiable, "[Error]: Unnecessary toValueGrid() call or toGradGrid() call");
    public:
        using ScalarType = T::PlainType;
    };

    template<Grid T>
    class Traits<GradGrid<T>> : public Traits<ValueGrid<T>> {};
}
