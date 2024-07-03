/*
 * Copyright 2023-2024 WeiBo He.
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
    template<class GridType>
    class RealGrid : public RValueGrid<RealGrid<GridType>> {
        using Base = RValueGrid<RealGrid<GridType>>;
        const GridType& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit RealGrid(const RValueGrid<GridType>& g_) : g(g_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).getReal(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class GridType>
    class ImagGrid : public RValueGrid<ImagGrid<GridType>> {
        using Base = RValueGrid<ImagGrid<GridType>>;
        const GridType& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit ImagGrid(const RValueGrid<GridType>& g_) : g(g_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).getImag(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class GridType>
    class NormGrid : public RValueGrid<NormGrid<GridType>> {
        using Base = RValueGrid<NormGrid<GridType>>;
        const GridType& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit NormGrid(const RValueGrid<GridType>& g_) : g(g_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).norm(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class GridType>
    class ValueGrid : public RValueGrid<ValueGrid<GridType>> {
        using Base = RValueGrid<ValueGrid<GridType>>;
        const GridType& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit ValueGrid(const RValueGrid<GridType>& g_) : g(g_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).getValue(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class GridType>
    class GradGrid : public RValueGrid<GradGrid<GridType>> {
        using Base = RValueGrid<GradGrid<GridType>>;
        const GridType& g;
    public:
        using typename Base::Index3D;
        using typename Base::ScalarType;
    public:
        explicit GradGrid(const RValueGrid<GridType>& g_) : g(g_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(Index3D index) const { return g.calc(index).getGrad(); }
        [[nodiscard]] size_t getDimX() const noexcept { return g.getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return g.getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return g.getDimZ(); }
    };

    template<class GridType>
    [[nodiscard]] inline RealGrid<GridType> toRealGrid(const RValueGrid<GridType>& v) {
        return RealGrid<GridType>{v};
    }

    template<class GridType>
    [[nodiscard]] inline ImagGrid<GridType> toImagGrid(const RValueGrid<GridType>& v) {
        return ImagGrid<GridType>{v};
    }

    template<class GridType>
    [[nodiscard]] inline NormGrid<GridType> toNormGrid(const RValueGrid<GridType>& v) {
        return NormGrid<GridType>{v};
    }

    template<class GridType>
    [[nodiscard]] inline ValueGrid<GridType> toValueGrid(const RValueGrid<GridType>& v) {
        return ValueGrid<GridType>{v};
    }

    template<class GridType>
    [[nodiscard]] inline GradGrid<GridType> toGradGrid(const RValueGrid<GridType>& v) {
        return GradGrid<GridType>{v};
    }
}

namespace Physica {
    using namespace Core;

    template<class GridType>
    class Traits<RealGrid<GridType>> {
    public:
        using ScalarType = typename GridType::ScalarType::RealType;
    };

    template<class GridType>
    class Traits<ImagGrid<GridType>> : public Traits<RealGrid<GridType>> {};

    template<class GridType>
    class Traits<NormGrid<GridType>> : public Traits<RealGrid<GridType>> {};

    template<class GridType>
    class Traits<ValueGrid<GridType>> {
        using T = typename GridType::ScalarType;
        static_assert(T::isDifferentiable, "[Error]: Unnecessary toValueGrid() call or toGradGrid() call");
    public:
        using ScalarType = typename T::PlainType;
    };

    template<class GridType>
    class Traits<GradGrid<GridType>> : public Traits<ValueGrid<GridType>> {};
}
