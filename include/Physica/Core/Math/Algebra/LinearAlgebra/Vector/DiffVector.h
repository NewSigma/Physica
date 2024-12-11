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

#include "Physica/Core/MultiPrecision/Diff.h"
#include "DenseVector.h"

namespace Physica::Core {
    template<Scalar T, int Order, size_t Length, class Allocator>
    class DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator> : public ContinuousVector<DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>> {
        using This = DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>;
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
        using ValueVector = DenseVector<T, Length>;
        using GradType = ScalarType::GradType;
        using GradVector = std::conditional<Order == 1, ValueVector, DenseVector<GradType, Length>>::type;

        ValueVector values;
        GradVector grads;
    public:
        DenseVector() = default;
        explicit DenseVector(size_t length);
        DenseVector(size_t length, const ScalarType& init);
        DenseVector(std::initializer_list<ScalarType> list);
        DenseVector(ValueVector values_, GradVector grads_) noexcept;
        template<Vector V>
        DenseVector(const V& v);
        DenseVector(const This&) = default;
        DenseVector(This&&) noexcept = default;
        ~DenseVector() = default;
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
        template<RandomGenerator R> inline void random_uniform();
        template<RandomGenerator R> inline void random_normal();
        template<class Distribution, RandomGenerator R>
        inline void random_any(Distribution& dist);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] size_t getLength() const noexcept { return values.getLength(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const noexcept;
        /* Static members */
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_uniform(size_t len);
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_normal(size_t len);
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist);
        [[nodiscard]] static auto linspace(T from, T to, size_t count);
    };
    ////////////////////////////////////////////////////////////////////////////////////
    template<Scalar T, int Order>
    class Diff<VectorND<T>, DiffMode::Reverse, Order>
            : public RValueVector<Diff<VectorND<T>, DiffMode::Reverse, Order>> {
        using VectorType = VectorND<T>;
        using This = Diff<VectorType, DiffMode::Reverse, Order>;
        using Base = RValueVector<This>;
        using TracerType = DiffTracer<T, Order>;
        using SegmentType = TraceSegment<T, Order>;
    public:
        using ScalarType = Base::ScalarType;
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
        template<RandomGenerator R> inline void random_uniform();
        template<RandomGenerator R> inline void random_normal();
        template<class Distribution, RandomGenerator R>
        inline void random_any(Distribution& dist);
        void swap(Diff& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return traceSeg.getLength(); }
        [[nodiscard]] const VectorType& getValue() const noexcept { return traceSeg.getValues(); }
        [[nodiscard]] const VectorType& getGrad() const noexcept { return traceSeg.getGrads(); }
        /* Static members */
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_uniform(size_t len);
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_normal(size_t len);
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist);
        /* Friends */
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, int Order, size_t Length, class Allocator>
    class Traits<Core::DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>> {
        static_assert(!T::isDiffable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, DiffMode::Forward, Order>;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<Scalar T, int Order>
    class Traits<Core::Diff<VectorND<T>, Core::DiffMode::Reverse, Order>> : public Traits<VectorND<T>> {
        static_assert(!T::isDiffable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Diff<T, DiffMode::Reverse, Order>;
    };
}

#include "DiffVectorImpl/DiffVectorImpl.h"
#include "DiffVectorImpl/Iterator.h"
