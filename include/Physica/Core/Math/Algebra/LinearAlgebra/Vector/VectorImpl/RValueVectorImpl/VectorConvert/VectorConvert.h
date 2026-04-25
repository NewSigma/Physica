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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

namespace Physica {
    template<class V>
    class ImagVector : public RValueVector<ImagVector<V>> {
        using This = ImagVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> v;
    public:
        explicit ImagVector(V&& v_) : v(std::forward<V>(v_)) {}
        ImagVector(const This&) = default;
        ImagVector(This&&) noexcept = default;
        ~ImagVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).imag(); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
    };

    template<class V>
    auto ImagVector<V>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v).values().imags();
    }

    template<class V>
    class NormVector : public RValueVector<NormVector<V>> {
        using This = NormVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> v;
    public:
        explicit NormVector(V&& v_) : v(std::forward<V>(v_)) {}
        NormVector(const This&) = default;
        NormVector(This&&) noexcept = default;
        ~NormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const { return v.calc(s).norm(); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
    };

    template<class V>
    auto NormVector<V>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v).values().norms();
    }

    template<class V>
    class ValueVector : public RValueVector<ValueVector<V>> {
        using This = ValueVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<V> v;
    public:
        explicit ValueVector(V&& v_) : v(std::forward<V>(v_)) {}
        ValueVector(const This&) = default;
        ValueVector(This&&) noexcept = default;
        ~ValueVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).value(); }
        [[nodiscard]] T calc_value(size_t s) const { return calc(s); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
    };

    template<class V, int GradOrder>
    class GradVector : public RValueVector<GradVector<V, GradOrder>> {
        using This = GradVector<V, GradOrder>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> v;
    public:
        explicit GradVector(V&& v_) : v(std::forward<V>(v_)) {}
        GradVector(const This&) = default;
        GradVector(This&&) noexcept = default;
        ~GradVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).template grad<GradOrder>(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return calc(s).value(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
    };

    template<class V, int MaskOrder>
    class GradMaskVector : public RValueVector<GradMaskVector<V, MaskOrder>> {
        static_assert(MaskOrder < std::remove_cvref_t<V>::ScalarType::Order, "[Error]: We should return ref to original vector instead of GradMaskVector, this is a bug");

        using This = GradMaskVector<V, MaskOrder>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> v;
    public:
        GradMaskVector(V&& v_) : v(std::forward<V>(v_)) {}
        GradMaskVector(const This&) = default;
        GradMaskVector(This&&) noexcept = default;
        ~GradMaskVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).template grad_mask<MaskOrder>(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<V>::getSizeAtCompile(); }
    };
}

namespace Physica {
    template<class V>
    class Traits<ImagVector<V>> : public Traits<RealVectorR<V>> {};

    template<class V>
    class Traits<NormVector<V>> : public Traits<RealVectorR<V>> {};

    template<class V>
    class Traits<ValueVector<V>> {
        using V1 = std::remove_cvref_t<V>;
    public:
        using ScalarType = V1::ScalarType::ValueType;
    };

    template<class V, int GradOrder>
    class Traits<GradVector<V, GradOrder>> {
        using V1 = std::remove_cvref_t<V>;
        static_assert(V1::ScalarType::isDiffable(), "[Error]: Redundant GradVector");
    public:
        using ScalarType = Internal::GradTypeHelper<typename V1::ScalarType, GradOrder>::Type;
    };

    template<class V, int MaskOrder>
    class Traits<GradMaskVector<V, MaskOrder>> {
        using V1 = std::remove_cvref_t<V>;
        using U = V1::ScalarType;
        using ValueType = typename U::ValueType;
        static_assert(U::isDiffable(), "[Error]: Redundant GradMaskVector");
    public:
        using ScalarType = std::conditional<MaskOrder == 0, ValueType, Diff<ValueType, U::Mode, MaskOrder>>::type;
    };
}
