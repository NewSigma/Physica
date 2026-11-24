/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/IndexVar/Index.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<class T, int... Dims>
    class ArrayND {
        using This = ArrayND<T, Dims...>;
        constexpr static bool StaticShape = sizeof...(Dims) != 1;
        static_assert(sizeof...(Dims) > 0, "[Error]: Dims is not specified");
    public:
        constexpr static int NDim = [](int x, auto...) consteval static noexcept -> int {
            return StaticShape ? sizeof...(Dims) : x;
        }(Dims...);

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

        using IndexType = Index<NDim>;
        using ArrayType = Helper<T>::Type;
        using ShapeType = std::conditional<StaticShape, Empty, IndexType>::type;

        ArrayType arr;
        [[no_unique_address]] ShapeType shape;
    public:
        ArrayND() = default;
        explicit ArrayND(IndexType shape_, auto&&... args);
        explicit ArrayND(std::same_as<size_t> auto... dims);
        ArrayND(const This&) = default;
        ArrayND(This&&) noexcept = default;
        ~ArrayND() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] auto& operator[](this auto&&, const IndexType& indices) noexcept;
        [[nodiscard]] auto& operator[](this auto&&, std::same_as<size_t> auto... dims) noexcept;
        [[nodiscard]] bool operator==(this const This& lhs, const This& rhs);
        [[nodiscard]] bool operator!=(this const This& lhs, const This& rhs) = default;
        /* Operations */
        void resize(IndexType shape_, auto&&... args);
        void resize(std::same_as<size_t> auto... dims);
        void reserve(size_t size) noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        [[nodiscard]] size_t toIndex1D(const IndexType& indices) const noexcept;
        [[nodiscard]] IndexType toIndexND(size_t index) const noexcept;
        void forND(this auto&& self, std::invocable<T&, IndexType> auto func);

        void zeros() noexcept;
        void junk() noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr IndexType getShape() const noexcept;
        [[nodiscard]] constexpr size_t dim(int index) const noexcept;
        [[nodiscard, gnu::returns_nonnull]] auto* data(this auto&&) noexcept;
        [[nodiscard, gnu::returns_nonnull]] auto* data_ptr(this auto&&, const IndexType& indices) noexcept;
        [[nodiscard]] auto&& asArray(this auto&&) noexcept;
        [[nodiscard]] constexpr size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] bool empty() const noexcept { return arr.empty(); }
        /* Static members */
        [[nodiscard]] constexpr static int ndim() noexcept { return NDim; }
        [[nodiscard]] consteval static bool isStaticShape() noexcept { return StaticShape; }
        [[nodiscard]] consteval static IndexType getShapeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ static size_t toSize(const IndexType& shape) noexcept;
        [[nodiscard]] static This zeros(IndexType shape);
        [[nodiscard]] static This zeros(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        [[nodiscard]] static This junk(IndexType shape);
        [[nodiscard]] static This junk(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        [[nodiscard]] static This read(IndexType shape, const T* __restrict p) noexcept;
        /* Friends */
        friend class device_obj<This>;
    };
}

#include "ArrayImpl/ArrayNDImpl.h"
