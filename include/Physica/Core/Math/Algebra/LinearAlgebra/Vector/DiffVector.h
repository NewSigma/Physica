/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Scalar/Diff.h" // IWYU pragma: export
#include "DenseVector.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    class DenseVector<Diff<T, Mode, Order>, Length, Allocator>
            : public CompactVector<DenseVector<Diff<T, Mode, Order>, Length, Allocator>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<DenseVector<Diff<T, Mode, Order>, Length, Allocator>>, Empty>::type {
        using This = DenseVector<Diff<T, Mode, Order>, Length, Allocator>;
        using Base = CompactVector<This>;
    public:
        using device_obj_type = device_obj<This>;
        using typename Base::ScalarType;
        template<Scalar U>
        using rebind_scalar = DenseVector<U, Length>;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    private:
        using ValueVector = DenseVector<T, Length>;
        using GradType = ScalarType::GradType;
        using GradVector = std::conditional<Order == 1, ValueVector, DenseVector<GradType, Length, Allocator>>::type;

        ValueVector v;
        GradVector g;
    public:
        DenseVector() = default;
        explicit DenseVector(size_t length) noexcept;
        DenseVector(size_t length, T init);
        DenseVector(size_t length, ScalarType init) requires(isForwardDiff());
        DenseVector(std::initializer_list<T> list);
        DenseVector(std::initializer_list<ScalarType> list) requires(isForwardDiff());
        DenseVector(ValueVector v_, GradVector g_) noexcept;
        explicit(isReverseDiff()) DenseVector(const Vector auto& src); // Force explicit conversion to avoid misuse
        DenseVector(const This&) = default;
        DenseVector(This&&) noexcept = default;
        ~DenseVector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        template<size_t I>
        [[nodiscard]] constexpr auto get(this auto&&) noexcept;
        [[nodiscard]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard]] constexpr auto end(this auto&&) noexcept;

        using Base::resize;
        void resize(this auto&, size_t size);

        template<RNG R> void random_uniform(this auto&);
        template<RNG R> void random_normal(this auto&);
        template<RNG R>
        void random_any(this auto&, auto& distribution);

        void zero_grad(this auto&);
        void junk(this auto&);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] auto data_handle(this auto&&) noexcept;

        [[nodiscard]] auto&& values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto&& grads(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t len, auto& distribution);
        [[nodiscard]] static This linspace(T from, T to, size_t count);
        [[nodiscard]] static This junk(size_t len);
        /* Friends */
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    class Traits<DenseVector<Diff<T, Mode, Order>, Length, Allocator>> : public Traits<VectorND<T>> {
        static_assert(!T::isDiffable(), "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, Mode, Order>;
    };
}

#include "DiffVectorImpl/DiffVectorImpl.h" // IWYU pragma: export
