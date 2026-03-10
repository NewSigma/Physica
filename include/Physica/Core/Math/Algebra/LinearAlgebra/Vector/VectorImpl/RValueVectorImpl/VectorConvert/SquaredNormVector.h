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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/Utils/Unroller.h"

namespace Physica {
    template<class V>
    class SquaredNormVector : public RValueVector<SquaredNormVector<V>> {
        using This = SquaredNormVector<V>;
        using Base = RValueVector<This>;
        using ComplexType = V::ScalarType;
    public:
        using Base::SizeAtCompile;
        using Base::isReverseDiff;
        using typename Base::PacketType;
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
        [[nodiscard]] Pack packet(size_t index) const noexcept;
        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index, size_t count) const noexcept;
        void reverse(const Scalar auto& grad) const noexcept;

        [[nodiscard]] CoDiff<T> sum() const noexcept;

        [[nodiscard]] auto values() const noexcept { return v.values().squaredNorms(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class V>
    template<Packet Pack>
    Pack SquaredNormVector<V>::packet(size_t index) const noexcept {
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

    template<class V>
    template<Packet Pack>
    Pack SquaredNormVector<V>::packet(size_t index, size_t count) const noexcept {
        assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
        if constexpr (isComplexV) {
            static_assert(!isReverseDiff, "[Error]: Not implemented");
            constexpr size_t MaxSize = BestPacket<ComplexType, Dynamic>::Size;
            constexpr size_t Size = Pack::size();
            if constexpr (Size <= MaxSize) {
                using PacketType = SIMD<ComplexType, Size>;
                const auto x2 = PacketType::asComplex(v.template packet<PacketType>(index, count).squaredNorm());
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
                    x2 = PacketType::asComplex(v.template packet<PacketType>(index, count).squaredNorm());

                if (count <= HalfSize)
                    return Pack(x2.real(), SIMD<T, HalfSize>(0));
                const auto y2 = PacketType::asComplex(v.template packet<PacketType>(index + HalfSize, count - HalfSize).squaredNorm());
                return Pack(x2.real(), y2.real());
            }
        }
        else
            return square(v.template packet<Pack>(index, count));
    }

    template<class V>
    void SquaredNormVector<V>::reverse(const Scalar auto& grad) const noexcept {
        static_assert(isReverseDiff);
        v.reverse(Tv(2) * grad * v.values());
    }

    template<class V>
    auto SquaredNormVector<V>::sum() const noexcept -> CoDiff<T> {
        if constexpr (!isReverseDiff && !isComplexV && Internal::EnableSIMD<V>::value) {
            assert(getLength() != 0 && "[Error]: Sum of a empty vector is not well defined");
            auto view = v.view();
            auto unroller = Unroller<PacketType, HostDevAttr::NumUnrollDefault>();
            size_t i = unroller.loop_recursive([ite = view.begin()](PacketType buffer, size_t index) noexcept {
                auto pack = (ite + index).template load<PacketType>();
                return fma(pack, pack, buffer);
            }, getLength());

            T result = unroller.sum().sum();
            for (auto ite = view.begin() + i; ite < view.end(); ++ite) {
                T x = *ite;
                result = fma(x, x, result);
            }
            return result;
        }
        else
            return Base::sum();
    }
}

namespace Physica {
    template<class V>
    class Traits<SquaredNormVector<V>> : public Traits<RealVectorR<V>> {};
}
