/*
 * Copyright 2021-2026 Weibo He.
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

#include "../CompactMatrix.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    class Transpose<M> : public CompactMatrix<Transpose<M>> {
        using This = Transpose<M>;
        using Base = CompactMatrix<This>;
    public:        
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
    public:
        Transpose(M&& mat_) : mat(std::forward<M>(mat_)) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const;
        void assign_mkl(Matrix auto&& target) const;

        [[nodiscard]] auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return mat.getCol(); }
        [[nodiscard]] size_t getCol() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    template<ExecutePolicy P>
    void Transpose<M>::assign(Matrix auto&& target) const {
        constexpr bool LargeMatrix = Base::getSizeAtCompile() == Dynamic;
        constexpr bool UseMKL = HasMKL() && Internal::EnableLAPACK<M, decltype(target)>::value && MatrixMajor::isSameMajor<M, decltype(target)>();
        if constexpr (LargeMatrix && UseMKL) {
            if (Base::getSize() <= 16)
                Base::template assign<P>(target);
            else
                assign_mkl(target);
        }
        else
            Base::template assign<P>(target);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    auto&& Transpose<M>::transpose(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    auto Transpose<M>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().values().transpose();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    auto Transpose<M>::data(this auto&& self) noexcept {
        return self.mat.data();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ consteval size_t Transpose<M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ consteval size_t Transpose<M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ consteval int Transpose<M>::getMajor() noexcept {
        constexpr static int OtherMajor = MatrixMajor::isColMatrix<M>() ? MatrixMajor::Row : MatrixMajor::Col;
        return MatrixMajor::isBothMajor<M>() ? MatrixMajor::BothMajor : OtherMajor;
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    class Transpose<V> : public CompactMatrix<Transpose<V>> {
        using This = Transpose<V>;
        using Base = CompactMatrix<This>;
    public:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<V> vec;
    public:
        explicit Transpose(V&& vec_) : vec(std::forward<V>(vec_)) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getCol() const noexcept { return vec.getLength(); }
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::Row; }
    };

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    auto&& Transpose<V>::transpose(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    auto Transpose<V>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().values().transpose();
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    auto Transpose<V>::data(this auto&& self) noexcept {
        return self.vec.data();
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ consteval size_t Transpose<V>::getRowAtCompile() noexcept {
        return 1;
    }

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ consteval size_t Transpose<V>::getColAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }
}

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    class Traits<Transpose<M>> {
    public:
        using ScalarType = std::remove_cvref_t<M>::ScalarType;
    };

    template<Vector V> requires(std::remove_cvref_t<V>::isCompact())
    class Traits<Transpose<V>> {
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
    };
}

#ifdef PHYSICA_MKL
    #include "Transpose_MKL.h"
#endif
