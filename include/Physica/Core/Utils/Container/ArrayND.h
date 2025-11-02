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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    using Index2D = Array<size_t, 2>;
    using Index3D = Array<size_t, 3>;
    using Index4D = Array<size_t, 4>;
    using Index5D = Array<size_t, 5>;
    using IndexND = Array<size_t>;

    template<class T, int Dim = Dynamic>
    class ArrayND {
        static_assert(Dim >= 0, "[Error]: Invalid Dim");
        using This = ArrayND<T, Dim>;
        using ShapeType = Array<size_t, Dim>;
        
        template<class U>
        struct Helper {
            using Type = Array<T>;
        };

        template<Scalar U>
        struct Helper<U> {
            using Type = DenseVector<T>;
        };
    private:
        using ArrayType = Helper<T>::Type;

        ArrayType arr;
        ShapeType shape;
    public:
        ArrayND() = default;
        explicit ArrayND(ShapeType shape_, auto&&... args);
        explicit ArrayND(size_t dim0, auto... dims);
        ArrayND(const This&) = default;
        ArrayND(This&&) noexcept = default;
        ~ArrayND() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T& operator()(const ShapeType& indices);
        [[nodiscard]] const T& operator()(const ShapeType& indices) const;
        [[nodiscard]] T& operator()(size_t dim0, auto... dims);
        [[nodiscard]] const T& operator()(size_t dim0, auto... dims) const;
        /* Operations */
        void resize(ShapeType shape_, auto&&... args);
        void resize(size_t dim0, auto... dims);
        void zeros() noexcept;

        [[nodiscard]] size_t toIndex1D(const ShapeType& indices) const noexcept;
        [[nodiscard]] ShapeType toIndexND(size_t index) const noexcept;

        void forND(std::invocable<T&, ShapeType> auto func);
        void forND(std::invocable<const T&, ShapeType> auto func) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto& asArray() noexcept { return arr; }
        [[nodiscard]] const auto& asArray() const noexcept { return arr; }
        [[nodiscard]] T* data_ptr(const ShapeType& indices) noexcept;
        [[nodiscard]] const T* data_ptr(const ShapeType& indices) const noexcept;
        [[nodiscard]] const ShapeType& getShape() const noexcept { return shape; }
        [[nodiscard]] size_t getShape(int dim) const noexcept { return shape[dim]; }
        [[nodiscard]] int getDim() const noexcept { return shape.getLength(); }
        [[nodiscard]] size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] bool empty() const noexcept { return arr.empty(); }
        /* Static members */
        [[nodiscard]] static size_t toSize(const ShapeType& shape) noexcept;
        [[nodiscard]] static size_t toIndex1D(const ShapeType& shape, const ShapeType& indices) noexcept;
        [[nodiscard]] static ShapeType toIndexND(const ShapeType& shape, size_t index) noexcept;
    private:
        static ShapeType toShape(auto... dims) noexcept;
        static void toShapeImpl(ShapeType& arr, int count, size_t dim0, auto... dims) noexcept;
        static void toShapeImpl(ShapeType&, int) noexcept {}
    };

    template<size_t Dim, class Functor>
    void forND(const Array<size_t, Dim>& shape, Functor func) {
        static_assert(Dim == 3, "[Error]: Not implemented");
        for (size_t x = 0; x < shape[0]; ++x)
            for (size_t y = 0; y < shape[1]; ++y)
                for (size_t z = 0; z < shape[2]; ++z)
                    func(Index3D{x, y, z});
    }
}

#include "ArrayImpl/ArrayNDImpl.h"
