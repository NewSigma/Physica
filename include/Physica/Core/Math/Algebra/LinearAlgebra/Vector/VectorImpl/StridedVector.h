/*
 * Copyright 2026 Weibo He.
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

#include "LValueVector.h"

namespace Physica {
    template<class Derived>
    class StridedVector : public LValueVector<Derived> {
        using Base = LValueVector<Derived>;
        using This = StridedVector<Derived>;

        template<Vector> class View;
    protected:
        using typename Base::T;
        using typename Base::Tr;
    public:
        ~StridedVector() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        Derived& operator=(Scalar auto x) noexcept;
        using Base::operator=;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] constexpr auto view(this auto&&) noexcept;

        [[nodiscard]] CoDiff<Tr> norm1() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm1_base() const noexcept;
        [[nodiscard]] Tr norm1_mkl() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm2() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm2_base() const noexcept;
        [[nodiscard]] Tr norm2_mkl() const noexcept;

        void zeros() noexcept;
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        /* Getters */
        [[nodiscard]] constexpr size_t getStride() const noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStrided() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastPacket() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept { return Dynamic; }
    protected:
        StridedVector() = default;
        StridedVector(const This&) = default;
        StridedVector(This&&) noexcept = default;
    };

    template<class Derived>
    constexpr auto StridedVector<Derived>::view(this auto&& self) noexcept {
        return View<std::remove_reference_t<decltype(self)>>(self);
    }

    template<class Derived>
    constexpr size_t StridedVector<Derived>::getStride() const noexcept {
        if constexpr (Derived::getStrideAtCompile() != Dynamic)
            return Derived::getStrideAtCompile();
        else
            return Base::getDerived().getStride();
    }

    template<class Derived>
    auto StridedVector<Derived>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedVector<Derived>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedVector<Derived>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength());
        return self.data_handle() + index * self.getStride();
    }

    template<class Derived>
    __host__ __device__ consteval bool StridedVector<Derived>::isCompact() noexcept {
        return Derived::getStrideAtCompile() == 1;
    }

    template<class Derived>
    __host__ __device__ consteval bool StridedVector<Derived>::isFastPacket() noexcept {
        return isCompact() || Derived::getStrideAtCompile() == 2;
    }
}

#include "StridedVectorImpl/StridedVectorImpl.h"
#include "StridedVectorImpl/View.h"
#ifdef PHYSICA_MKL
    #include "StridedVectorImpl/StridedVector_MKL.h"
#endif
