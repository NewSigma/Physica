/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Scalar/Diff.h" // IWYU pragma: export
#include "DenseVector.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    class DenseVector<Diff<T, Mode, Order>, Length, Allocator>
            : public ContinuousVector<DenseVector<Diff<T, Mode, Order>, Length, Allocator>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<DenseVector<Diff<T, Mode, Order>, Length, Allocator>>, PlainStruct<void>>::type {
        using This = DenseVector<Diff<T, Mode, Order>, Length, Allocator>;
        using Base = ContinuousVector<This>;
    public:
        using device_obj_type = device_obj<This>;
        using typename Base::ScalarType;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
        using IterF = PtrIteratorF<This>;
        using IterR = PtrIteratorR<This>;
        using IterCF = PtrIteratorF<const This>;
        using IterCR = PtrIteratorR<const This>;
    private:
        using ValueVector = DenseVector<T, Length>;
        using GradType = ScalarType::GradType;
        using GradVector = std::conditional<Order == 1, ValueVector, DenseVector<GradType, Length, Allocator>>::type;

        ValueVector v;
        GradVector g;
    public:
        DenseVector() = default;
        explicit DenseVector(size_t length);
        DenseVector(size_t length, T init);
        DenseVector(size_t length, ScalarType init) requires(isForwardDiff);
        DenseVector(std::initializer_list<T> list);
        DenseVector(std::initializer_list<ScalarType> list) requires(isForwardDiff);
        DenseVector(ValueVector v_, GradVector g_) noexcept;
        explicit(isReverseDiff) DenseVector(const Vector auto& src); // Force explicit conversion to avoid misuse
        DenseVector(const This&) = default;
        DenseVector(This&&) noexcept = default;
        ~DenseVector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const;
        using Base::operator=;
        using Base::operator[];
        /* Iterators */
        using Base::begin;
        using Base::cbegin;
        using Base::end;
        using Base::cend;
        using Base::rbegin;
        using Base::crbegin;
        using Base::rend;
        using Base::crend;
        /* Operations */
        void zero_grad();

        using Base::resize;
        void resize(size_t size);

        template<RNG R = Random<>> void random_uniform();
        template<RNG R = Random<>> void random_normal();
        template<RNG R = Random<>>
        void random_any(auto& distribution);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] ConstPtrTy data_ptr(size_t index) const noexcept;

        [[nodiscard]] const auto& values() const noexcept { return v; }
        [[nodiscard]] auto& values() noexcept { return v; }
        [[nodiscard]] const auto& grads() const noexcept { return g; }
        [[nodiscard]] auto& grads() noexcept { return g; }
        /* Static members */
        template<RNG R = Random<>>
        [[nodiscard]] static This random_uniform(size_t len);
        template<RNG R = Random<>>
        [[nodiscard]] static This random_normal(size_t len);
        template<RNG R = Random<>>
        [[nodiscard]] static This random_any(size_t len, auto& distribution);
        [[nodiscard]] static auto linspace(T from, T to, size_t count);
        /* Friends */
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    class Traits<DenseVector<Diff<T, Mode, Order>, Length, Allocator>> : public Traits<VectorND<T>> {
        static_assert(!T::isDiffable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, Mode, Order>;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}

#include "DiffVectorImpl/DiffVectorImpl.h" // IWYU pragma: export
#include "DiffVectorImpl/Iterator.h" // IWYU pragma: export
