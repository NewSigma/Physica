/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Scalar/Diff.h" // IWYU pragma: export
#include "DenseTensor.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    class DenseTensor<Diff<T, Mode, Order>, Dims...>
            : public CompactTensor<DenseTensor<Diff<T, Mode, Order>, Dims...>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<DenseTensor<Diff<T, Mode, Order>, Dims...>>, Empty>::type {
        using This = DenseTensor<Diff<T, Mode, Order>, Dims...>;
        using Base = CompactTensor<This>;
    public:
        using typename Base::IndexType;
        using typename Base::ScalarType;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
        template<Scalar U>
        using rebind_scalar = DenseTensor<U, Dims...>;
    private:
        using ValueTensor = DenseTensor<T, Dims...>;
        using GradType = ScalarType::GradType;
        using GradTensorType = std::conditional<Order == 1, ValueTensor, DenseTensor<GradType, Dims...>>::type;

        ValueTensor v;
        GradTensorType g;
    public:
        DenseTensor() = default;
        explicit DenseTensor(IndexType shape);
        DenseTensor(IndexType shape, T init);
        DenseTensor(IndexType shape, ScalarType init) requires(isForwardDiff());
        DenseTensor(ValueTensor v_, GradTensorType g_);
        explicit(isReverseDiff()) DenseTensor(const Tensor auto& x);
        DenseTensor(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        DenseTensor(const This&) = default;
        DenseTensor(This&&) noexcept = default;
        ~DenseTensor() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(this auto&, IndexType shape);
        void reserve(this auto&, size_t size) noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        template<RNG R>
        void random_uniform(this auto&);
        template<RNG R>
        void random_normal(this auto&);
        template<RNG R>
        void random_any(this auto&, auto& distribution);

        void zero_grad(this auto&);
        void junk(this auto&);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto data_handle(this auto&&) noexcept;
        [[nodiscard]] constexpr IndexType getShape() const noexcept { return v.getShape(); }
        [[nodiscard]] constexpr size_t dim(int index) const noexcept { return v.dim(index); }
        [[nodiscard]] constexpr size_t getSize() const noexcept { return v.getSize(); }
        [[nodiscard]] bool empty() const noexcept { return v.empty(); }

        [[nodiscard]] auto&& values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto&& grads(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticShape() noexcept { return ValueTensor::isStaticShape(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return ValueTensor::getSizeAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static auto getShapeAtCompile() noexcept -> IndexType { return ValueTensor::getShapeAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static auto getStrideAtCompile() noexcept -> IndexType { return ValueTensor::getStrideAtCompile(); }

        [[nodiscard]] static This zeros(IndexType shape);
        [[nodiscard]] static This zeros(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        [[nodiscard]] static This junk(IndexType shape);
        [[nodiscard]] static This junk(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        template<RNG R>
        [[nodiscard]] static This random_uniform(IndexType shape);
        template<RNG R>
        [[nodiscard]] static This random_uniform(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        template<RNG R>
        [[nodiscard]] static This random_normal(IndexType shape);
        template<RNG R>
        [[nodiscard]] static This random_normal(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0);
        template<RNG R>
        [[nodiscard]] static This random_any(IndexType shape, auto& distribution);
        /* Friends */
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    class Traits<DenseTensor<Diff<T, Mode, Order>, Dims...>> : public Traits<DenseTensor<T, Dims...>> {
        static_assert(!T::isDiffable(), "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, Mode, Order>;
    };
}

#include "TensorImpl/DiffDenseTensorImpl.h"
