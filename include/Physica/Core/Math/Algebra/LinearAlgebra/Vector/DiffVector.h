/*
 * Copyright 2024 Weibo He.
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

#include <Physica/Core/MultiPrecision/Diff.h>
#include "Vector.h"

namespace Physica::Core {
    template<class T, int Order, size_t Length, class Allocator>
    class Vector<Diff<T, DiffMode::Forward, Order>, Length, Allocator> : public ContinuousVector<Vector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>> {
        using This = Vector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
        using FIteType = FIterator<This>;
        using RIteType = RIterator<This>;
        using CFIteType = FIterator<const This>;
        using CRIteType = RIterator<const This>;
    private:
        using ValueVector = Vector<T, Length>;
        using GradType = typename ScalarType::GradType;
        using GradVector = typename std::conditional<Order == 1, ValueVector, Vector<GradType, Length>>::type;

        ValueVector values;
        GradVector grads;
    public:
        Vector() = default;
        explicit Vector(size_t length);
        Vector(size_t length, const ScalarType& init);
        Vector(std::initializer_list<ScalarType> list);
        Vector(ValueVector values_, GradVector grads_) noexcept;
        template<class Derived>
        Vector(const RValueVector<Derived>& v);
        Vector(const This&) = default;
        Vector(This&&) noexcept = default;
        ~Vector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Iterators */
        __host__ __device__ FIteType begin() noexcept { return FIteType(data()); }
        __host__ __device__ CFIteType begin() const noexcept { return cbegin(); }
        __host__ __device__ CFIteType cbegin() const noexcept { return CFIteType(data()); }
        __host__ __device__ FIteType end() noexcept { return FIteType(data() + getLength()); }
        __host__ __device__ CFIteType end() const noexcept { return cend(); }
        __host__ __device__ CFIteType cend() const noexcept { return CFIteType(data() + getLength()); }
        __host__ __device__ RIteType rbegin() noexcept { return RIteType(data() + getLength() - 1); }
        __host__ __device__ CRIteType rbegin() const noexcept { return crbegin(); }
        __host__ __device__ CRIteType crbegin() const noexcept { return CRIteType(data() + getLength() - 1); }
        __host__ __device__ RIteType rend() noexcept { return RIteType(data() - 1); }
        __host__ __device__ CRIteType rend() const noexcept { return crend(); }
        __host__ __device__ CRIteType crend() const noexcept { return CRIteType(data() - 1); }
        /* Operations */
        inline void resize(size_t size);
        template<class RandomType> inline void random_uniform(RandomType& gen);
        template<class RandomType> inline void random_normal(RandomType& gen);
        template<class Distribution, class RandomType>
        inline void random_any(Distribution& dist, RandomType& gen);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] size_t getLength() const noexcept { return values.getLength(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const noexcept;
        /* Static members */
        template<class RandomType>
        [[nodiscard]] inline static This random_uniform(size_t len, RandomType& gen);
        template<class RandomType>
        [[nodiscard]] inline static This random_normal(size_t len, RandomType& gen);
        template<class Distribution, class RandomType>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist, RandomType& gen);
        [[nodiscard]] static auto linspace(ScalarType from, ScalarType to, size_t count);
    };
    ////////////////////////////////////////////////////////////////////////////////////
    template<class T, int Order>
    class Diff<Vector<T>, DiffMode::Reverse, Order>
            : public RValueVector<Diff<Vector<T>, DiffMode::Reverse, Order>> {
        using VectorType = Vector<T>;
        using This = Diff<VectorType, DiffMode::Reverse, Order>;
        using Base = RValueVector<This>;
        using TracerType = DiffTracer<T, Order>;
        using SegmentType = TraceSegment<T, Order>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Diff();
        Diff(size_t length);
        Diff(VectorType values);
        Diff(const Diff&) = default;
        Diff(Diff&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        Diff& operator=(const Diff&) = default;
        Diff& operator=(Diff&&) noexcept = default;
        /* Operations */
        template<class RandomGenerator> inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator> inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(Diff& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return traceSeg.getLength(); }
        [[nodiscard]] const VectorType& getValue() const noexcept { return traceSeg.getValues(); }
        [[nodiscard]] const VectorType& getGrad() const noexcept { return traceSeg.getGrads(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist, RandomGenerator& gen);
        /* Friends */
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<class T, int Order, size_t Length, class Allocator>
    class Traits<Core::Vector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, DiffMode::Forward, Order>;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class T, int Order>
    class Traits<Core::Diff<Vector<T>, Core::DiffMode::Reverse, Order>> : public Traits<Vector<T>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, DiffMode::Reverse, Order>;
    };
}

#include "DiffVectorImpl/DiffVectorImpl.h"
#include "DiffVectorImpl/Iterator.h"
