/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "TensorImpl/CompactTensor.h"

namespace Physica {
    template<Scalar T, int... Dims>
    class DenseTensor : public CompactTensor<DenseTensor<T, Dims...>> {
        using This = DenseTensor<T, Dims...>;
        using Base = CompactTensor<This>;
    public:
        using Base::NDim;
        using typename Base::IndexType;
    private:
        ArrayND<T, Dims...> storage;
    public:
        DenseTensor() = default;
        explicit DenseTensor(ArrayND<T, Dims...> storage) noexcept;
        DenseTensor(IndexType shape, auto&&... args);
        DenseTensor(std::same_as<size_t> auto... dims);
        DenseTensor(const Tensor auto& x);
        DenseTensor(const This&) = default;
        DenseTensor(This&&) noexcept = default;
        ~DenseTensor() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using Base::resize;
        void resize(this auto&, IndexType shape);
        void reserve(this auto&, size_t size) noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        using Base::random_normal;
        using Base::random_uniform;
        using Base::random_any;
        void zeros(this auto&) noexcept;
        void junk(this auto&) noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto data_handle(this auto&&) noexcept;
        [[nodiscard]] constexpr IndexType getShape() const noexcept;
        [[nodiscard]] constexpr size_t dim(int index) const noexcept;
        [[nodiscard]] auto&& asArray(this auto&&) noexcept;
        [[nodiscard]] constexpr size_t getSize() const noexcept;
        /* Static members */
        using Base::getStrideAtCompile;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticShape() noexcept { return ArrayND<T, Dims...>::isStaticShape(); }
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static auto getShapeAtCompile() noexcept -> IndexType { return ArrayND<T, Dims...>::getShapeAtCompile(); }
        [[nodiscard]] static This zeros(IndexType shape);
        [[nodiscard]] static This zeros(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        [[nodiscard]] static This junk(IndexType shape);
        [[nodiscard]] static This junk(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        template<RNG R>
        static DenseTensor random_uniform(IndexType shape);
        template<RNG R>
        static DenseTensor random_uniform(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        template<RNG R>
        static DenseTensor random_normal(IndexType shape);
        template<RNG R>
        static DenseTensor random_normal(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        template<RNG R>
        [[nodiscard]] static This random_any(IndexType shape, auto& distribution);
        [[nodiscard]] static This generate(std::invocable<IndexType> auto fn, IndexType shape);
        [[nodiscard]] static This generate(std::invocable<IndexType> auto fn, std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        [[nodiscard]] static This read(IndexType shape, const T* __restrict p) noexcept;
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T, int... Dims>
    void swap(DenseTensor<T, Dims...>& x, DenseTensor<T, Dims...>& y) noexcept {
        x.swap(y);
    }

    template<Scalar T> using Tensor3D = DenseTensor<T, 3>;
    template<Scalar T> using Tensor4D = DenseTensor<T, 4>;
    template<Scalar T> using TensorND = DenseTensor<T>;
}

namespace Physica {
    template<Scalar T, int... Dims>
    class Traits<DenseTensor<T, Dims...>> {
        static_assert(sizeof...(Dims) > 0, "[Error]: Dims is not specified");
        static_assert(!Diffable<T>, "[Error]: Use diffable tensor instead");
    public:
        using ScalarType = T;
        constexpr static int NDim = ArrayND<T, Dims...>::NDim;
    };
}

#include "TensorImpl/DenseTensorImpl.h"
