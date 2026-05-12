/*
 * Copyright 2024-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix M, Vector V>
    class GEMV : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M&& mat, V&& vec);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Vector auto&& target) const noexcept;
        void assign_mkl(Vector auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto&& target) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        template<int Size>
        [[nodiscard]] auto packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] auto packet(size_t index, size_t count) const noexcept;

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        [[nodiscard]] auto grads(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ constexpr static bool isFastPacket() noexcept;
    };

    template<Matrix M, Vector V>
    GEMV<M, V>::GEMV(M&& mat, V&& vec) : mat(std::forward<M>(mat)), vec(std::forward<V>(vec)) {}

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto&& target) const noexcept {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        using T1 = M1::ScalarType;
        using T2 = V1::ScalarType;
        if constexpr (std::same_as<T1, T2> && M1::isCompact() && HasMKL() && Internal::EnableLAPACK<V1, decltype(target)>::value) {
            constexpr size_t Size = std::max(M1::getColAtCompile(), V1::getSizeAtCompile());
            constexpr size_t Threhold = 16;
            if constexpr (Size == Dynamic) {
                if (vec.getLength() > Threhold)
                    assign_mkl(target);
                else
                    assign_base<P>(target);
            }
            else if constexpr (Size > Threhold)
                assign_mkl(target);
            else
                assign_base<P>(target);
        }
        else
            assign_base<P>(target);
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_base(Vector auto&& target) const noexcept {
        target.assert_assign(*this);
        if constexpr (std::remove_cvref_t<V>::isFastAssign()) {
            DenseVector<T, Base::getSizeAtCompile(target)> buffer = getRHS();
            (getLHS() * std::move(buffer)).template assign<P>(target);
        }
        else if constexpr (MatrixMajor::isColMatrix<M>()) {
            size_t length = vec.getLength();
            (mat.col(0) * vec.calc(0)).template assign<P>(target);
            for (size_t i = 1; i < length; ++i)
                (mat.col(i) * vec.calc(i)).template assign_add<P>(target);
        }
        else {
            size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                target[i] = calc(i);
        }
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_add(Vector auto&& target) const noexcept {
        if constexpr (MatrixMajor::isColMatrix<M>()) {
            size_t length = vec.getLength();
            for (size_t i = 0; i < length; ++i)
                (mat.col(i) * vec.calc(i)).template assign_add<P>(target);
        }
        else {
            size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                target[i] += calc(i);
        }
    }

    template<Matrix M, Vector V>
    auto GEMV<M, V>::calc(size_t index) const -> CoDiff<T> {
        return mat.row(index) * vec;
    }

    template<Matrix M, Vector V>
    template<int Size>
    auto GEMV<M, V>::packet(size_t index) const noexcept {
        if constexpr (isFastPacket() && !Base::isDiffable()) {
            size_t length = vec.getLength();
            auto unroller = Unroller<SIMD<T, Size>, HostDevAttr::NumUnrollGEMV>();
            [[maybe_unused]] size_t i = unroller.loop_recursive([this, index](SIMD<T, Size> buffer, size_t iSize) noexcept {
                size_t i = iSize / Size;
                return fma(mat.col(i).template packet<Size>(index), SIMD<T, Size>(vec.calc(i)), buffer);
            }, length * Size);
            assert(i == length * Size && "[Error]: Unroller is buggy");
            return unroller.sum();
        }
        else
            return Base::template packet<Size>(index);
    }

    template<Matrix M, Vector V>
    template<int Size>
    auto GEMV<M, V>::packet(size_t index, size_t count) const noexcept {
        if constexpr (isFastPacket() && !Base::isDiffable()) {
            size_t length = vec.getLength();
            auto unroller = Unroller<SIMD<T, Size>, HostDevAttr::NumUnrollGEMV>();
            [[maybe_unused]] size_t i = unroller.loop_recursive([this, index, count](SIMD<T, Size> buffer, size_t iSize) noexcept {
                size_t i = iSize / Size;
                return fma(mat.col(i).template packet<Size>(index, count), SIMD<T, Size>(vec.calc(i)), buffer);
            }, length * Size);
            assert(i == length * Size && "[Error]: Unroller is buggy");
            return unroller.sum();
        }
        else
            return Base::template packet<Size>(index, count);
    }

    template<Matrix M, Vector V>
    void GEMV<M, V>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        Base::assert_assign(grad);
        const auto& g = grad.values();
        if constexpr (ReverseDiff<M>) {
            if constexpr (MatrixMajor::isRowMatrix<M>()) {
                for (size_t i = 0; i < mat.getRow(); ++i)
                    mat.row(i).reverse(g[i] * vec.values());
            }
            else {
                for (size_t i = 0; i < mat.getCol(); ++i)
                    mat.col(i).reverse(g * vec.calc_value(i));
            }
        }

        if constexpr (ReverseDiff<V>) {
            if constexpr (MatrixMajor::isRowMatrix<M>()) {
                for (size_t i = 0; i < getLength(); ++i)
                    vec.reverse(mat.values().row(i) * g.calc(i));
            }
            else
                vec.reverse(mat.values().transpose() * g);
        }
    }

    template<Matrix M, Vector V>
    auto GEMV<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V>
    auto GEMV<M, V>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (Base::isForwardDiff()) {
            if constexpr (ForwardDiff<M> && ForwardDiff<V>)
                return std::forward<Self>(self).getLHS().grads() * std::forward<Self>(self).getRHS().values()
                     + std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().grads();
            else if constexpr (ForwardDiff<M>)
                return std::forward<Self>(self).getLHS().grads() * std::forward<Self>(self).getRHS().values();
            else
                return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().grads();
        }
        else
            return self.Base::grads();
    }

    template<Matrix M, Vector V>
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V>
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Matrix M, Vector V>
    __host__ __device__ consteval size_t GEMV<M, V>::getSizeAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<Matrix M, Vector V>
    __host__ __device__ consteval bool GEMV<M, V>::isFastAssign() noexcept {
        return MatrixMajor::isColMatrix<M>();
    }

    template<Matrix M, Vector V>
    __host__ __device__ constexpr bool GEMV<M, V>::isFastPacket() noexcept {
        return MatrixMajor::isColMatrix<M>();
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        static_assert(M1::getColAtCompile() == V1::getSizeAtCompile() || M1::getColAtCompile() == Dynamic || V1::getSizeAtCompile() == Dynamic,
                "Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename V1::ScalarType>::Type;
    };
}

#ifdef PHYSICA_MKL
    #include "GEMV_MKL.h"
#endif
