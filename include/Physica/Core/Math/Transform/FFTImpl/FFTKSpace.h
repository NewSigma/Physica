/*
 * Copyright 2023 WeiBo He.
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
    template<class Derived> class FFTKSpace;

    namespace Internal {
        template<class Derived>
        class Traits<FFTKSpace<Derived>> {
            using T = typename Traits<Derived>::ScalarType;
        public:
            using ScalarType = typename T::ComplexType;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class Derived>
    class FFTKSpace
            : public Utils::CRTPBase<Derived, 1>
            , public ContinuousVector<FFTKSpace<Derived>> {
        using This = FFTKSpace<Derived>;
        using Base = Utils::CRTPBase<Derived, 1>;
        using VectorBase = ContinuousVector<FFTKSpace<Derived>>;
    public:
        using ScalarType = typename Internal::Traits<This>::ScalarType;
    private:
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
    public:
        /* Operators */
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator[](size_t index);
        [[nodiscard]] const ScalarType& operator[](size_t index) const;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Derived::rSizeToKSize(Base::getDerived().getRSpaceSize()); }
        [[nodiscard]] size_t getSize() const noexcept { return getLength(); }
    };

    template<class Derived>
    inline typename FFTKSpace<Derived>::ScalarType& FFTKSpace<Derived>::operator[](size_t index) {
        return Base::getDerived().complex_buffer[index];
    }

    template<class Derived>
    inline const typename FFTKSpace<Derived>::ScalarType& FFTKSpace<Derived>::operator[](size_t index) const {
        return Base::getDerived().complex_buffer[index];
    }
}
