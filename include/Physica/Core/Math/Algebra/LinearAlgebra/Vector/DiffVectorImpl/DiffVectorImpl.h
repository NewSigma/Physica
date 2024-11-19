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

namespace Physica::Core {
    template<Scalar T, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::DenseVector(size_t length) : values(length), grads(length) {}

    template<Scalar T, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::DenseVector(size_t length, const ScalarType& init)
            : values(length, init.getValue()), grads(length, init.getGrad()) {}

    template<Scalar T, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::DenseVector(std::initializer_list<ScalarType> list) : DenseVector(list.size()) {
        size_t i = 0;
        for (auto& elem : list) {
            values[i] = elem.getValue();
            grads[i] = elem.getGrad();
            i += 1;
        }
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::DenseVector(ValueVector values_, GradVector grads_) noexcept
            : values(std::move(values_)), grads(std::move(grads_)) {}

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<Vector V>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::DenseVector(const V& v) : DenseVector(v.getLength()) {
        v.assignTo(*this);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    inline void DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::resize(size_t size) {
        values.resize(size);
        grads.resize(size);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<class RandomType>
    inline void DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::random_uniform(RandomType& gen) {
        values.random_uniform(gen);
        grads = ScalarType(0);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<class RandomType>
    inline void DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::random_normal(RandomType& gen) {
        values.random_normal(gen);
        grads = ScalarType(0);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<class Distribution, class RandomType>
    inline void DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::random_any(Distribution& dist, RandomType& gen) {
        values.random_any(dist, gen);
        grads = ScalarType(0);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    void DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    __host__ __device__ inline typename DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::PtrTy
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::data_ptr(size_t index) noexcept {
        assert(index < getLength() && "[Error]: Index out of range");
        return PtrTy(values.data_ptr(index), grads.data_ptr(index));
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    __host__ __device__ inline typename DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::ConstPtrTy
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::data_ptr(size_t index) const noexcept {
        return const_cast<This&>(*this).data_ptr(index);
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<class RandomType>
    inline DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::random_uniform(size_t len, RandomType& gen) {
        This result(len);
        result.random_uniform(gen);
        return result;
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<class RandomType>
    inline DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::random_normal(size_t len, RandomType& gen) {
        This result(len);
        result.random_normal(gen);
        return result;
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    template<class Distribution, class RandomType>
    inline DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>
    DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::random_any(size_t len, Distribution& dist, RandomType& gen) {
        This result(len);
        result.random_any(dist, gen);
        return result;
    }

    template<Scalar T, int Order, size_t Length, class Allocator>
    auto DenseVector<Diff<T, DiffMode::Forward, Order>, Length, Allocator>::linspace(ScalarType from, ScalarType to, size_t count) {
        auto values = ValueVector::linspace(from.getValue(), to.getValue(), count);
        auto grads = GradVector::linspace(from.getGrad(), to.getGrad(), count);
        return This(std::move(values), std::move(grads));
    }
    ////////////////////////////////////////////////////////////////////////////////////
    template<Scalar T, int Order>
    Diff<VectorND<T>, DiffMode::Reverse, Order>::Diff()
            : traceSeg(TracerType::getInstance().pushSegment()) {}

    template<Scalar T, int Order>
    Diff<VectorND<T>, DiffMode::Reverse, Order>::Diff(size_t length)
            : traceSeg(TracerType::getInstance().pushSegment(length)) {}

    template<Scalar T, int Order>
    Diff<VectorND<T>, DiffMode::Reverse, Order>::Diff(VectorType values)
            : traceSeg(TracerType::getInstance().pushSegment(std::move(values))) {}

    template<Scalar T, int Order>
    template<class RandomGenerator>
    inline void Diff<VectorND<T>, DiffMode::Reverse, Order>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getLength(), gen);
    }

    template<Scalar T, int Order>
    template<class RandomGenerator>
    inline void Diff<VectorND<T>, DiffMode::Reverse, Order>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getLength(), gen);
    }

    template<Scalar T, int Order>
    template<class Distribution, class RandomGenerator>
    inline void Diff<VectorND<T>, DiffMode::Reverse, Order>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getLength(), dist, gen);
    }

    template<Scalar T, int Order>
    inline typename Diff<VectorND<T>, DiffMode::Reverse, Order>::ScalarType
    Diff<VectorND<T>, DiffMode::Reverse, Order>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        return traceSeg[index];
    }

    template<Scalar T, int Order>
    template<class RandomGenerator>
    inline Diff<VectorND<T>, DiffMode::Reverse, Order>
    Diff<VectorND<T>, DiffMode::Reverse, Order>::random_uniform(size_t len, RandomGenerator& gen) {
        return This(VectorType::random_uniform(len, gen));
    }

    template<Scalar T, int Order>
    template<class RandomGenerator>
    inline Diff<VectorND<T>, DiffMode::Reverse, Order>
    Diff<VectorND<T>, DiffMode::Reverse, Order>::random_normal(size_t len, RandomGenerator& gen) {
        return This(VectorType::random_normal(len, gen));
    }

    template<Scalar T, int Order>
    template<class Distribution, class RandomGenerator>
    inline Diff<VectorND<T>, DiffMode::Reverse, Order>
    Diff<VectorND<T>, DiffMode::Reverse, Order>::random_any(size_t len, Distribution& dist, RandomGenerator& gen) {
        return This(VectorType::random_any(len, dist, gen));
    }
}
