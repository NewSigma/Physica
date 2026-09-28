/*
 * Copyright 2023-2026 Weibo He.
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

#include "../RValueTensor.h"

namespace Physica {
    template<class X>
    class RealTensor : public RValueTensor<RealTensor<X>> {
        using This = RealTensor<X>;
        using Base = RValueTensor<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
    public:
        explicit RealTensor(X&& tensor_) : tensor(std::forward<X>(tensor_)) {}
        RealTensor(const This&) = default;
        RealTensor(This&&) noexcept = default;
        ~RealTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(const IndexType& index) const { return tensor.calc(index).real(); }

        [[nodiscard]] decltype(auto) values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return tensor.getShape(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getSizeAtCompile(); }
    };

    template<class X>
    decltype(auto) RealTensor<X>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), X>(self.tensor).values().reals();
    }

    template<class X>
    class ImagTensor : public RValueTensor<ImagTensor<X>> {
        using This = ImagTensor<X>;
        using Base = RValueTensor<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
    public:
        explicit ImagTensor(X&& tensor_) : tensor(std::forward<X>(tensor_)) {}
        ImagTensor(const This&) = default;
        ImagTensor(This&&) noexcept = default;
        ~ImagTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(const IndexType& index) const { return tensor.calc(index).imag(); }

        [[nodiscard]] decltype(auto) values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return tensor.getShape(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getSizeAtCompile(); }
    };

    template<class X>
    decltype(auto) ImagTensor<X>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), X>(self.tensor).values().imags();
    }

    template<class X>
    class SquaredNormTensor : public RValueTensor<SquaredNormTensor<X>> {
        using This = SquaredNormTensor<X>;
        using Base = RValueTensor<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
    public:
        explicit SquaredNormTensor(X&& tensor_) : tensor(std::forward<X>(tensor_)) {}
        SquaredNormTensor(const This&) = default;
        SquaredNormTensor(This&&) noexcept = default;
        ~SquaredNormTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(const IndexType& index) const { return tensor.calc(index).squaredNorm(); }

        [[nodiscard]] decltype(auto) values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return tensor.getShape(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getSizeAtCompile(); }
    };

    template<class X>
    decltype(auto) SquaredNormTensor<X>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), X>(self.tensor).values().squaredNorms();
    }

    template<class X>
    class NormTensor : public RValueTensor<NormTensor<X>> {
        using This = NormTensor<X>;
        using Base = RValueTensor<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
    public:
        explicit NormTensor(X&& tensor_) : tensor(std::forward<X>(tensor_)) {}
        NormTensor(const This&) = default;
        NormTensor(This&&) noexcept = default;
        ~NormTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(const IndexType& index) const { return tensor.calc(index).norm(); }

        [[nodiscard]] decltype(auto) values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return tensor.getShape(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getSizeAtCompile(); }
    };

    template<class X>
    decltype(auto) NormTensor<X>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), X>(self.tensor).values().norms();
    }

    template<class X>
    class ValueTensor : public RValueTensor<ValueTensor<X>> {
        using This = ValueTensor<X>;
        using Base = RValueTensor<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
    public:
        explicit ValueTensor(X&& tensor_) : tensor(std::forward<X>(tensor_)) {}
        ValueTensor(const This&) = default;
        ValueTensor(This&&) noexcept = default;
        ~ValueTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(const IndexType& index) const { return tensor.calc(index).value(); }
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return tensor.getShape(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getSizeAtCompile(); }
    };

    template<class X, int GradOrder>
    class GradTensor : public RValueTensor<GradTensor<X, GradOrder>> {
        using This = GradTensor<X, GradOrder>;
        using Base = RValueTensor<This>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<X> tensor;
    public:
        explicit GradTensor(X&& tensor_) : tensor(std::forward<X>(tensor_)) {}
        GradTensor(const This&) = default;
        GradTensor(This&&) noexcept = default;
        ~GradTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) calc(const IndexType& index) const { return tensor.calc(index).template grad<GradOrder>(); }
        /* Getters */
        [[nodiscard]] IndexType getShape() const noexcept { return tensor.getShape(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<X>::getSizeAtCompile(); }
    };
}

namespace Physica {
    template<class X>
    class Traits<RealTensor<X>> {
        using X1 = std::remove_cvref_t<X>;
    public:
        using ScalarType = X1::ScalarType::RealType;
        constexpr static int NDim = X1::NDim;
    };

    template<class X>
    class Traits<ImagTensor<X>> : public Traits<RealTensor<X>> {};

    template<class X>
    class Traits<SquaredNormTensor<X>> : public Traits<RealTensor<X>> {};

    template<class X>
    class Traits<NormTensor<X>> : public Traits<RealTensor<X>> {};

    template<class X>
    class Traits<ValueTensor<X>> {
        using X1 = std::remove_cvref_t<X>;
    public:
        using ScalarType = X1::ScalarType::ValueType;
        constexpr static int NDim = X1::NDim;
    };

    template<class X, int GradOrder>
    class Traits<GradTensor<X, GradOrder>> {
        using X1 = std::remove_cvref_t<X>;
        static_assert(X1::isDiffable(), "[Error]: Redundant GradTensor");
    public:
        using ScalarType = X1::ScalarType::template GradWithOrder<GradOrder>;
        constexpr static int NDim = X1::NDim;
    };
}
