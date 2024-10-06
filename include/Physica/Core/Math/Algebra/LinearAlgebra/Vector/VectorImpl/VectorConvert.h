/*
 * Copyright 2023-2024 Weibo He.
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
    template<class VectorType>
    class RealVector : public RValueVector<RealVector<VectorType>> {
        using This = RealVector<VectorType>;
        using Base = RValueVector<This>;

        const VectorType& v;
    public:
        explicit RealVector(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        RealVector(const This&) = delete;
        RealVector(This&&) noexcept = delete;
        ~RealVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        typename Base::ScalarType calc(size_t s) const { return v.calc(s).getReal(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class ImagVector : public RValueVector<ImagVector<VectorType>> {
        using This = ImagVector<VectorType>;
        using Base = RValueVector<This>;

        const VectorType& v;
    public:
        explicit ImagVector(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        ImagVector(const This&) = delete;
        ImagVector(This&&) noexcept = delete;
        ~ImagVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        typename Base::ScalarType calc(size_t s) const { return v.calc(s).getImag(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class SquaredNormVector : public RValueVector<SquaredNormVector<VectorType>> {
        using This = SquaredNormVector<VectorType>;
        using Base = RValueVector<This>;

        const VectorType& v;
    public:
        explicit SquaredNormVector(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        SquaredNormVector(const This&) = delete;
        SquaredNormVector(This&&) noexcept = delete;
        ~SquaredNormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        typename Base::ScalarType calc(size_t s) const { return v.calc(s).squaredNorm(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class NormVector : public RValueVector<NormVector<VectorType>> {
        using This = NormVector<VectorType>;
        using Base = RValueVector<This>;

        const VectorType& v;
    public:
        explicit NormVector(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        NormVector(const This&) = delete;
        NormVector(This&&) noexcept = delete;
        ~NormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        typename Base::ScalarType calc(size_t s) const { return v.calc(s).norm(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class ValueVector : public RValueVector<ValueVector<VectorType>> {
        using This = ValueVector<VectorType>;
        using Base = RValueVector<This>;

        const VectorType& v;
    public:
        explicit ValueVector(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        ValueVector(const This&) = delete;
        ValueVector(This&&) noexcept = delete;
        ~ValueVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        typename Base::ScalarType calc(size_t s) const { return v.calc(s).getValue(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class GradVector : public RValueVector<GradVector<VectorType>> {
        using This = GradVector<VectorType>;
        using Base = RValueVector<This>;

        const VectorType& v;
    public:
        explicit GradVector(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        GradVector(const This&) = delete;
        GradVector(This&&) noexcept = delete;
        ~GradVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        typename Base::ScalarType calc(size_t s) const { return v.calc(s).getGrad(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    [[nodiscard]] inline auto toRealVector(const RValueVector<VectorType>& v) {
        return RealVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toImagVector(const RValueVector<VectorType>& v) {
        return ImagVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toSquaredNormVector(const RValueVector<VectorType>& v) {
        return SquaredNormVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toNormVector(const RValueVector<VectorType>& v) {
        return NormVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toValueVector(const RValueVector<VectorType>& v) {
        return ValueVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toGradVector(const RValueVector<VectorType>& v) {
        return GradVector<VectorType>{v};
    }
}

namespace Physica {
    using namespace Core;

    template<class VectorType>
    class Traits<RealVector<VectorType>> {
    public:
        using ScalarType = typename VectorType::ScalarType::RealType;
        constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = Traits<VectorType>::FastPacket;
    };

    template<class VectorType>
    class Traits<ImagVector<VectorType>> : public Traits<RealVector<VectorType>> {};

    template<class VectorType>
    class Traits<SquaredNormVector<VectorType>> : public Traits<RealVector<VectorType>> {};

    template<class VectorType>
    class Traits<NormVector<VectorType>> : public Traits<RealVector<VectorType>> {};

    template<class VectorType>
    class Traits<ValueVector<VectorType>> {
        using T = typename VectorType::ScalarType;
        static_assert(T::isDifferentiable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = typename T::PlainScalar;
        constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<class VectorType>
    class Traits<GradVector<VectorType>> : public Traits<ValueVector<VectorType>> {};
}
