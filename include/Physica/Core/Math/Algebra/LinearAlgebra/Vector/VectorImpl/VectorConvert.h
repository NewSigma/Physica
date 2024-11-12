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
    public:
        using typename Base::ScalarType;
    private:
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
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).real(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class ImagVector : public RValueVector<ImagVector<VectorType>> {
        using This = ImagVector<VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
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
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).imag(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class SquaredNormVector : public RValueVector<SquaredNormVector<VectorType>> {
        using This = SquaredNormVector<VectorType>;
        using Base = RValueVector<This>;
        using ComplexType = typename VectorType::ScalarType;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    private:
        constexpr static bool isComplexV = VectorType::isComplex;

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
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).squaredNorm(); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            if constexpr (isComplexV) {
                static_assert(!isReverseDiff, "[Error]: Not implemented");
                constexpr size_t MaxSize = BestPacket<ComplexType, Dynamic>::Size;
                constexpr size_t Size = AnyPacket::size();
                if constexpr (Size <= MaxSize) {
                    using PacketType = SIMD<ComplexType, Size>;
                    const auto reim = PacketType::asComplex(square(v.template packet<PacketType>(index).asReal()));
                    const auto imre = PacketType::asComplex(reim.swapRealImag());
                    return AnyPacket((reim + imre).real());
                }
                else {
                    constexpr size_t Size1 = Size / 2;
                    using PacketType = SIMD<ComplexType, Size1>;
                    const auto reim1 = PacketType::asComplex(square(v.template packet<PacketType>(index).asReal()));
                    const auto imre1 = PacketType::asComplex(reim1.swapRealImag());
                    const auto reim2 = PacketType::asComplex(square(v.template packet<PacketType>(index + Size1).asReal()));
                    const auto imre2 = PacketType::asComplex(reim2.swapRealImag());
                    return AnyPacket((reim1 + imre1).real(), (reim2 + imre2).real());
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
                    const auto reim = PacketType::asComplex(square(v.template packetPartial<PacketType>(index, count).asReal()));
                    const auto imre = PacketType::asComplex(reim.swapRealImag());
                    return AnyPacket((reim + imre).real());
                }
                else {
                    constexpr size_t Size1 = Size / 2;
                    using PacketType = SIMD<ComplexType, Size1>;
                    const bool flag = count <= Size1;
                    const size_t count1 = flag ? count : Size1;
                    const auto reim1 = PacketType::asComplex(square(v.template packetPartial<PacketType>(index, count1).asReal()));
                    const auto imre1 = PacketType::asComplex(reim1.swapRealImag());

                    if (flag)
                        return AnyPacket((reim1 + imre1).real(), SIMD<ScalarType, Size1>(0));
                    const size_t count2 = count - count1;
                    const auto reim2 = PacketType::asComplex(square(v.template packetPartial<PacketType>(index + Size / 2, count2).asReal()));
                    const auto imre2 = PacketType::asComplex(reim2.swapRealImag());
                    return AnyPacket((reim1 + imre1).real(), (reim2 + imre2).real());
                }
            }
            else
                return square(v.template packetPartial<AnyPacket>(index, count));
        }

        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class NormVector : public RValueVector<NormVector<VectorType>> {
        using This = NormVector<VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
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
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).norm(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class ValueVector : public RValueVector<ValueVector<VectorType>> {
        using This = ValueVector<VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
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
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).getValue(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType, int GradOrder>
    class GradVector : public RValueVector<GradVector<VectorType, GradOrder>> {
        using This = GradVector<VectorType, GradOrder>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
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
        [[nodiscard]] ScalarType calc(size_t s) const { return v.calc(s).template getGrad<GradOrder>(); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    [[nodiscard]] inline auto toRealVector(const RValueVector<VectorType>& v) noexcept {
        return RealVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toImagVector(const RValueVector<VectorType>& v) noexcept {
        return ImagVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toSquaredNormVector(const RValueVector<VectorType>& v) noexcept {
        return SquaredNormVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toNormVector(const RValueVector<VectorType>& v) noexcept {
        return NormVector<VectorType>{v};
    }

    template<class VectorType>
    [[nodiscard]] inline auto toValueVector(const RValueVector<VectorType>& v) noexcept {
        return ValueVector<VectorType>{v};
    }

    template<class VectorType, int GradOrder = 1>
    [[nodiscard]] inline auto toGradVector(const RValueVector<VectorType>& v) noexcept {
        return GradVector<VectorType, GradOrder>{v};
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
        using ScalarType = typename T::ValueType;
        constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<class VectorType, int GradOrder>
    class Traits<GradVector<VectorType, GradOrder>> {
        using T = typename VectorType::ScalarType;
        static_assert(T::isDifferentiable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = typename T::template GradRtnTy<GradOrder>;
        constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
