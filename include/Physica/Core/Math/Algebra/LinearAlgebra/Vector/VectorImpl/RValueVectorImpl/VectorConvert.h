/*
 * Copyright 2023-2025 Weibo He.
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

#include "../RValueVector.h"

namespace Physica {
    template<class V>
    class RealVector : public RValueVector<RealVector<V>> {
        using This = RealVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& v;
    public:
        explicit RealVector(const V& v_) : v(v_) {}
        RealVector(const This&) = default;
        RealVector(This&&) noexcept = default;
        ~RealVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).real(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s).real(); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V>
    class ImagVector : public RValueVector<ImagVector<V>> {
        using This = ImagVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& v;
    public:
        explicit ImagVector(const V& v_) : v(v_) {}
        ImagVector(const This&) = default;
        ImagVector(This&&) noexcept = default;
        ~ImagVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).imag(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s).imag(); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V>
    class SquaredNormVector : public RValueVector<SquaredNormVector<V>> {
        using This = SquaredNormVector<V>;
        using Base = RValueVector<This>;
        using ComplexType = V::ScalarType;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        constexpr static bool isComplexV = V::isComplex;

        const V& v;
    public:
        explicit SquaredNormVector(const V& v_) : v(v_) {}
        SquaredNormVector(const This&) = default;
        SquaredNormVector(This&&) noexcept = default;
        ~SquaredNormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] CoDiff<T> calc(size_t s) const { return v.calc(s).squaredNorm(); }

        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s).squaredNorm(); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            if constexpr (isComplexV) {
                static_assert(!isReverseDiff, "[Error]: Not implemented");
                constexpr size_t MaxSize = BestPacket<ComplexType, Dynamic>::Size;
                constexpr size_t Size = Pack::size();
                if constexpr (Size <= MaxSize) {
                    using PacketType = SIMD<ComplexType, Size>;
                    const auto x2 = PacketType::asComplex(v.template packet<PacketType>(index).squaredNorm());
                    return Pack(x2.real());
                }
                else {
                    constexpr size_t Size1 = Size / 2;
                    using PacketType = SIMD<ComplexType, Size1>;
                    const auto x2 = PacketType::asComplex(v.template packet<PacketType>(index).squaredNorm());
                    const auto y2 = PacketType::asComplex(v.template packet<PacketType>(index + Size1).squaredNorm());
                    return Pack(x2.real(), y2.real());
                }
            }
            else
                return square(v.template packet<Pack>(index));
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
            if constexpr (isComplexV) {
                static_assert(!isReverseDiff, "[Error]: Not implemented");
                constexpr size_t MaxSize = BestPacket<ComplexType, Dynamic>::Size;
                constexpr size_t Size = Pack::size();
                if constexpr (Size <= MaxSize) {
                    using PacketType = SIMD<ComplexType, Size>;
                    const auto x2 = PacketType::asComplex(v.template packetPartial<PacketType>(index, count).squaredNorm());
                    return Pack(x2.real());
                }
                else {
                    // We cannot finish the work in one run, separate the results into a low half and a high half.
                    constexpr size_t HalfSize = Size / 2;
                    using PacketType = SIMD<ComplexType, HalfSize>;
                    PacketType x2{};
                    if (count >= HalfSize)
                        x2 = PacketType::asComplex(v.template packet<PacketType>(index).squaredNorm());
                    else
                        x2 = PacketType::asComplex(v.template packetPartial<PacketType>(index, count).squaredNorm());

                    if (count <= HalfSize)
                        return Pack(x2.real(), SIMD<T, HalfSize>(0));
                    const auto y2 = PacketType::asComplex(v.template packetPartial<PacketType>(index + HalfSize, count - HalfSize).squaredNorm());
                    return Pack(x2.real(), y2.real());
                }
            }
            else
                return square(v.template packetPartial<Pack>(index, count));
        }

        void reverse(const Scalar auto& grad) const noexcept {
            static_assert(isReverseDiff);
            v.reverse(Tv(2) * grad * v.values());
        }

        [[nodiscard]] auto values() const noexcept { return v.values().squaredNorms(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V>
    class NormVector : public RValueVector<NormVector<V>> {
        using This = NormVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& v;
    public:
        explicit NormVector(const V& v_) : v(v_) {}
        NormVector(const This&) = default;
        NormVector(This&&) noexcept = default;
        ~NormVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] CoDiff<T> calc(size_t s) const { return v.calc(s).norm(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s).norm(); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V>
    class ValueVector : public RValueVector<ValueVector<V>> {
        using This = ValueVector<V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
    private:
        const V& v;
    public:
        explicit ValueVector(const V& v_) : v(v_) {}
        ValueVector(const This&) = default;
        ValueVector(This&&) noexcept = default;
        ~ValueVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] T calc(size_t s) const { return v.calc_value(s); }
        [[nodiscard]] T calc_value(size_t s) const { return calc(s); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V, int GradOrder>
    class GradVector : public RValueVector<GradVector<V, GradOrder>> {
        using This = GradVector<V, GradOrder>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& v;
    public:
        explicit GradVector(const V& v_) : v(v_) {}
        GradVector(const This&) = default;
        GradVector(This&&) noexcept = default;
        ~GradVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).template grad<GradOrder>(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return calc(s).value(); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V, int MaskOrder>
    class DiffMaskVector : public RValueVector<DiffMaskVector<V, MaskOrder>> {
        static_assert(MaskOrder < V::ScalarType::Order, "[Error]: We should return ref to original vector instead of DiffMaskVector, this is a bug");

        using This = DiffMaskVector<V, MaskOrder>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& v;
    public:
        DiffMaskVector(const V& v_) : v(v_) {}
        DiffMaskVector(const This&) = default;
        DiffMaskVector(This&&) noexcept = default;
        ~DiffMaskVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] T calc(size_t s) const { return v.calc(s).template mask<MaskOrder>(); }
        [[nodiscard]] Tv calc_value(size_t s) const { return v.calc_value(s); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<Vector V, int MaskOrder>
    [[nodiscard]] auto toDiffMaskVector(const RValueVector<V>& v) noexcept -> std::conditional<(MaskOrder < V::ScalarType::Order), DiffMaskVector<V, MaskOrder>, const V&>::type{
        return v.getDerived();
    }
}

namespace Physica {
    template<class V>
    class Traits<RealVector<V>> {
    public:
        using ScalarType = V::ScalarType::RealType;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = Traits<V>::FastPacket;
    };

    template<class V>
    class Traits<ImagVector<V>> : public Traits<RealVector<V>> {};

    template<class V>
    class Traits<SquaredNormVector<V>> : public Traits<RealVector<V>> {};

    template<class V>
    class Traits<NormVector<V>> : public Traits<RealVector<V>> {};

    template<class V>
    class Traits<ValueVector<V>> {
    public:
        using ScalarType = V::ScalarType::ValueType;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<class V, int GradOrder>
    class Traits<GradVector<V, GradOrder>> {
        static_assert(V::ScalarType::isDiffable, "[Error]: Unnecessary toGradVector() call");
    public:
        using ScalarType = Internal::GradTypeHelper<typename V::ScalarType, GradOrder>::Type;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<class V, int MaskOrder>
    class Traits<DiffMaskVector<V, MaskOrder>> {
        using U = V::ScalarType;
        using ValueType = typename U::ValueType;
        static_assert(U::isDiffable, "[Error]: Unnecessary toDiffMaskVector() call");
    public:
        using ScalarType = std::conditional<MaskOrder == 0, ValueType, Diff<ValueType, U::Mode, MaskOrder>>::type;
        constexpr static size_t SizeAtCompile = V::SizeAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
