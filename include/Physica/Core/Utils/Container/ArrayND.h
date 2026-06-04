/*
 * Copyright 2023-2026 Weibo He.
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
    template<class T, int... Dims>
    class ArrayND {
        using This = ArrayND<T, Dims...>;
        constexpr static bool StaticShape = sizeof...(Dims) != 1;
        static_assert(sizeof...(Dims) > 0, "[Error]: Dims is not specified");
    public:
        constexpr static int NDim = [](int x, auto...) consteval static noexcept -> int { return StaticShape ? sizeof...(Dims) : x; }(Dims...);
        constexpr static size_t SizeAtCompile = StaticShape ? (Dims * ...) : Dynamic;
        static_assert(NDim > 2, "[Error]: Invalid Dim");
    private:
        template<class U>
        struct Helper {
            using Type = Array<T, SizeAtCompile>;
        };

        template<Scalar U>
        struct Helper<U> {
            using Type = DenseVector<T, SizeAtCompile>;
        };

        using IndexType = Array<size_t, NDim>;
        using ArrayType = Helper<T>::Type;
        using ShapeType = std::conditional<StaticShape, Empty, IndexType>::type;

        ArrayType arr;
        [[no_unique_address]] ShapeType shape;
    public:
        ArrayND() = default;
        explicit ArrayND(IndexType shape_, auto&&... args);
        explicit ArrayND(std::integral auto... dims);
        ArrayND(const This&) = default;
        ArrayND(This&&) noexcept = default;
        ~ArrayND() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] auto& operator[](this auto&&, const IndexType& indices) noexcept;
        [[nodiscard]] auto& operator[](this auto&&, std::integral auto... dims) noexcept;
        /* Operations */
        void resize(IndexType shape_, auto&&... args);
        void resize(std::integral auto... dims);
        void zeros() noexcept;

        [[nodiscard]] size_t toIndex1D(const IndexType& indices) const noexcept;
        [[nodiscard]] IndexType toIndexND(size_t index) const noexcept;

        void forND(std::invocable<T&, IndexType> auto func);
        void forND(std::invocable<const T&, IndexType> auto func) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr int ndim() const noexcept { return NDim; }
        [[nodiscard]] constexpr size_t dim(int index) const noexcept;
        [[nodiscard]] IndexType getShape() const noexcept;
        [[nodiscard, gnu::returns_nonnull]] auto* data(this auto&&) noexcept;
        [[nodiscard, gnu::returns_nonnull]] auto* data_ptr(this auto&&, const IndexType& indices) noexcept;
        [[nodiscard]] auto& asArray() noexcept { return arr; }
        [[nodiscard]] const auto& asArray() const noexcept { return arr; }
        [[nodiscard]] size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] bool empty() const noexcept { return arr.empty(); }
        /* Static members */
        [[nodiscard]] static size_t toSize(const IndexType& shape) noexcept;
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
