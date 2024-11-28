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
    template<Vector T>
    class RealVector : public RValueVector<RealVector<T>> {
        using This = RealVector<T>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& v;
    public:
        explicit RealVector(const T& v_) : v(v_) {}
        RealVector(const This&) = delete;
        RealVector(This&&) noexcept = delete;
        ~RealVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).real(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<Vector T>
    class ImagVector : public RValueVector<ImagVector<T>> {
        using This = ImagVector<T>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& v;
    public:
        explicit ImagVector(const T& v_) : v(v_) {}
        ImagVector(const This&) = delete;
        ImagVector(This&&) noexcept = delete;
        ~ImagVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).imag(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<Vector T>
    class SquaredNormVector : public RValueVector<SquaredNormVector<T>> {
        using This = SquaredNormVector<T>;
        using Base = RValueVector<This>;
        using ComplexType = T::ScalarType;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    private:
        constexpr static bool isComplexV = T::isComplex;

        const T& v;
    public:
        explicit SquaredNormVector(const T& v_) : v(v_) {}
        SquaredNormVector(const This&) = delete;
        SquaredNormVector(This&&) noexcept = delete;
        ~SquaredNormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).squaredNorm(); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            if constexpr (isComplexV) {
                static_assert(!isReverseDiff, "[Error]: Not implemented");
                constexpr size_t MaxSize = BestPacket<ComplexType, Dynamic>::Size;
                constexpr size_t Size = AnyPacket::size();
                if constexpr (Size <= MaxSize) {
                    using PacketType = SIMD<ComplexType, Size>;
                    const auto x2 = PacketType::asComplex(v.template packet<PacketType>(index).squaredNorm());
                    return AnyPacket(x2.real());
                }
                else {
                    constexpr size_t Size1 = Size / 2;
                    using PacketType = SIMD<ComplexType, Size1>;
                    const auto x2 = PacketType::asComplex(v.template packet<PacketType>(index).squaredNorm());
                    const auto y2 = PacketType::asComplex(v.template packet<PacketType>(index + Size1).squaredNorm());
                    return AnyPacket(x2.real(), y2.real());
                }
            }
            else
                return square(v.template packet<AnyPacket>(index));
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            if constexpr (isComplexV) {
                static_assert(!isReverseDiff, "[Error]: Not implemented");
                constexpr size_t MaxSize = BestPacket<ComplexType, Dynamic>::Size;
                constexpr size_t Size = AnyPacket::size();
                if constexpr (Size <= MaxSize) {
                    using PacketType = SIMD<ComplexType, Size>;
                    const auto x2 = PacketType::asComplex(v.template packetPartial<PacketType>(index, count).squaredNorm());
                    return AnyPacket(x2.real());
                }
                else {
                    constexpr size_t Size1 = Size / 2;
                    using PacketType = SIMD<ComplexType, Size1>;
                    const bool flag = count <= Size1;
                    const size_t count1 = flag ? count : Size1;
                    const auto x2 = PacketType::asComplex(v.template packetPartial<PacketType>(index, count1).squaredNorm());

                    if (flag)
                        return AnyPacket(x2.real(), SIMD<ScalarType, Size1>(0));
                    const size_t count2 = count - count1;
                    const auto y2 = PacketType::asComplex(v.template packetPartial<PacketType>(index + Size / 2, count2).squaredNorm());
                    return AnyPacket(x2.real(), y2.real());
                }
            }
            else
                return square(v.template packetPartial<AnyPacket>(index, count));
        }

        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<Vector T>
    class NormVector : public RValueVector<NormVector<T>> {
        using This = NormVector<T>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& v;
    public:
        explicit NormVector(const T& v_) : v(v_) {}
        NormVector(const This&) = delete;
        NormVector(This&&) noexcept = delete;
        ~NormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).norm(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<Vector T>
    class ValueVector : public RValueVector<ValueVector<T>> {
        using This = ValueVector<T>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& v;
    public:
        explicit ValueVector(const T& v_) : v(v_) {}
        ValueVector(const This&) = delete;
        ValueVector(This&&) noexcept = delete;
        ~ValueVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).getValue(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<Vector T, int GradOrder>
    class GradVector : public RValueVector<GradVector<T, GradOrder>> {
        using This = GradVector<T, GradOrder>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& v;
    public:
        explicit GradVector(const T& v_) : v(v_) {}
        GradVector(const This&) = delete;
        GradVector(This&&) noexcept = delete;
        ~GradVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).template getGrad<GradOrder>(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<Vector T>
    [[nodiscard]] inline auto toRealVector(const T& v) noexcept {
        return RealVector<T>{v};
    }

    template<Vector T>
    [[nodiscard]] inline auto toImagVector(const T& v) noexcept {
        return ImagVector<T>{v};
    }

    template<Vector T>
    [[nodiscard]] inline auto toSquaredNormVector(const T& v) noexcept {
        return SquaredNormVector<T>{v};
    }

    template<Vector T>
    [[nodiscard]] inline auto toNormVector(const T& v) noexcept {
        return NormVector<T>{v};
    }

    template<Vector T>
    [[nodiscard]] inline auto toValueVector(const T& v) noexcept {
        return ValueVector<T>{v};
    }

    template<Vector T, int GradOrder = 1>
    [[nodiscard]] inline auto toGradVector(const T& v) noexcept {
        return GradVector<T, GradOrder>{v};
    }
}

namespace Physica {
    template<Vector T>
    class Traits<RealVector<T>> {
    public:
        using ScalarType = T::ScalarType::RealType;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = Traits<T>::FastPacket;
    };

    template<Vector T>
    class Traits<ImagVector<T>> : public Traits<RealVector<T>> {};

    template<Vector T>
    class Traits<SquaredNormVector<T>> : public Traits<RealVector<T>> {};

    template<Vector T>
    class Traits<NormVector<T>> : public Traits<RealVector<T>> {};

    template<Vector T>
    class Traits<ValueVector<T>> {
        static_assert(T::ScalarType::isDifferentiable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = T::ValueType;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<Vector T, int GradOrder>
    class Traits<GradVector<T, GradOrder>> {
        static_assert(T::ScalarType::isDifferentiable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = T::ScalarType::template GradRtnTy<GradOrder>;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
