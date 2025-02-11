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

namespace Physica {
    /**
     * FIterator(Forward)
     */
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    class FIterator<DenseVector<Diff<T, Mode, Order>, Length, Allocator>> {
        using ElemType = Diff<T, Mode, Order>;
        using Container = DenseVector<ElemType, Length, Allocator>;
        using This = FIterator<Container>;
        constexpr static bool isConst = std::is_const<Container>::value;
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = std::conditional<isConst, const ElemType, ElemType>::type;
        using pointer = ElemType::PtrTy;
        using reference = ElemType::RefTy;
        using iterator_category = std::random_access_iterator_tag;
    private:
        pointer p;
    public:
        __host__ __device__ explicit FIterator(pointer p) : p(p) {}
        __host__ __device__ FIterator(const This& ite) : p(ite.p) {}
        ~FIterator() = default;
        /* Operators */
        This& operator=(const This&) = default;
        [[nodiscard]] __host__ __device__ This operator+(difference_type n) const { return This(p + n); }
        [[nodiscard]] __host__ __device__ This operator-(difference_type n) const { return This(p - n); }
        [[nodiscard]] __host__ __device__ bool operator<(const This& ite) const { return p < ite.p; }
        [[nodiscard]] __host__ __device__ bool operator>(const This& ite) const { return p > ite.p; }
        [[nodiscard]] __host__ __device__ bool operator<=(const This& ite) const { return !(p > ite.p); }
        [[nodiscard]] __host__ __device__ bool operator>=(const This& ite) const { return !(p < ite.p); }
        [[nodiscard]] __host__ __device__ difference_type operator-(const This& ite) const { return p - ite.p; }
        [[nodiscard]] __host__ __device__ bool operator==(const This& ite) const noexcept { return p == ite.p; }
        [[nodiscard]] __host__ __device__ bool operator!=(const This& ite) const noexcept { return p != ite.p; }
        __host__ __device__ This& operator++() { ++p; return *this; }
        __host__ __device__ const This operator++(int) { return This(p++); }
        __host__ __device__ This& operator--() { --p; return *this; }
        [[nodiscard]] __host__ __device__ reference operator*() const { return *p; }
    };
    /**
     * RIterator(Reverse)
     */
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    class RIterator<DenseVector<Diff<T, Mode, Order>, Length, Allocator>> {
        using ElemType = Diff<T, Mode, Order>;
        using Container = DenseVector<ElemType, Length, Allocator>;
        using This = RIterator<Container>;
        constexpr static bool isConst = std::is_const<Container>::value;
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = std::conditional<isConst, const ElemType, ElemType>::type;
        using pointer = ElemType::PtrTy;
        using reference = ElemType::RefTy;
        using iterator_category = std::random_access_iterator_tag;
    private:
        pointer p;
    public:
        __host__ __device__ explicit RIterator(pointer p) : p(p) {}
        __host__ __device__ RIterator(const This& ite) : p(ite.p) {}
        ~RIterator() = default;
        /* Operators */
        This& operator=(const This&) = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& ite) const noexcept { return p == ite.p; }
        [[nodiscard]] __host__ __device__ bool operator!=(const This& ite) const noexcept { return p != ite.p; }
        __host__ __device__ This& operator++() { --p; return *this; }
        __host__ __device__ const This operator++(int)  { return This(p--); }
        [[nodiscard]] __host__ __device__ reference operator*() const { return *p; }
    };
}
