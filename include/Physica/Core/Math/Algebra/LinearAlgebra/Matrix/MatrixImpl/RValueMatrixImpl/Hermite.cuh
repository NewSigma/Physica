/*
 * Copyright 2025-2026 Weibo He.
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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<class T>
    struct remove_hermite<device_obj<Hermite<T>>> {
        using Type = device_obj<T>;
    };

    template<Matrix M>
    class device_obj<Hermite<M>> : public device_obj<RValueMatrix<Hermite<M>>> {
        using host_obj = Hermite<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        explicit device_obj(Ref mat_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] __host__ __device__ auto&& hermite(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return hermite().getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return hermite().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept { return hermite().getOrder(); }
    };

    template<Matrix M>
    device_obj<Hermite<M>>::device_obj(Ref mat_) : mat(asStruct(mat_)) {}

    template<Matrix M>
    __device__ auto device_obj<Hermite<M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        return hermite().calc(col, row, block).conjugate();
    }

    template<Matrix M>
    __host__ __device__ auto&& device_obj<Hermite<M>>::hermite(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<Hermite<M>>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).hermite().values();
    }

    template<Vector V>
    class device_obj<Hermite<V>> : public device_obj<RValueMatrix<Hermite<V>>> {
        using host_obj = Hermite<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<V>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> vec;
    public:
        explicit device_obj(Ref vec_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __device__ void assign(Matrix auto& target) const;

        using Base::calc;
        [[nodiscard]] __device__ T calc([[maybe_unused]] size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] __host__ __device__ auto&& hermite(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return hermite().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
    };

    template<Vector V>
    device_obj<Hermite<V>>::device_obj(Ref vec_) : vec(asStruct(vec_)) {}

    template<Vector V>
    __device__ void device_obj<Hermite<V>>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < hermite().getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(target.rowFromMajorMinor(0, i), target.colFromMajorMinor(0, i));
    }

    template<Vector V>
    __device__ auto device_obj<Hermite<V>>::calc([[maybe_unused]] size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        assert(row == 0);
        return hermite().calc(col, block).conjugate();
    }

    template<Vector V>
    __host__ __device__ auto&& device_obj<Hermite<V>>::hermite(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.vec.getDerived());
    }

    template<Vector V>
    __host__ __device__ auto device_obj<Hermite<V>>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).hermite().values();
    }

    template<Vector V>
    __host__ __device__ size_t device_obj<Hermite<V>>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return 1;
    }
}

namespace Physica {
    template <Matrix M>
    class Traits<device_obj<Hermite<M>>> : public Traits<Hermite<M>> {};

    template <Vector V>
    class Traits<device_obj<Hermite<V>>> : public Traits<Hermite<V>> {};
}

