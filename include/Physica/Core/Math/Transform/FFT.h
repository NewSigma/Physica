/*
 * Copyright 2020-2022 WeiBo He.
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

#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType, size_t Dim> class FFTImpl;
    }

    template<class ScalarType, size_t Dim = 1> class FFT;

    template<class ScalarType>
    class FFT<ScalarType, 1> {
        using Impl = Internal::FFTImpl<ScalarType, 1>;
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
        static constexpr bool isComplex = ScalarType::isComplex;
    private:
        Impl impl;
    public:
        FFT(size_t size_, const RealType& deltaT);
        FFT(const Vector<ScalarType>& data, const RealType& deltaT);
        FFT(const FFT&) = default;
        FFT(FFT&&) noexcept = default;
        ~FFT() = default;
        /* Operators */
        FFT& operator=(FFT fft) noexcept;
        ComplexType operator()(size_t i) { return impl(i); }
        /* Operations */
        void transform(const Vector<ScalarType>& data) { impl.transform(data); }
        void invTransform(const Vector<ComplexType>& data) { impl.invTransform(data); }
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return impl.getSize(); }
        [[nodiscard]] const RealType& getDeltaT() const noexcept { return impl.getDeltaT(); }
        [[nodiscard]] ComplexType getComponent(ssize_t index) const { return impl.getComponent(index); }
        [[nodiscard]] Vector<ComplexType> getComponents() const { return impl.getComponents(); }
        [[nodiscard]] RealType getDeltaFreq() const noexcept { return reciprocal(impl.getDeltaT() * impl.getSize()); }
        [[nodiscard]] ComplexType getFreqIntense(const RealType& freq) const noexcept;
        /* Helpers */
        void swap(FFT& fft) noexcept { impl.swap(fft.impl); }
    };

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(size_t size, const RealType& deltaT) : impl(size, deltaT) {}

    template<class ScalarType>
    FFT<ScalarType, 1>::FFT(const Vector<ScalarType>& data, const RealType& deltaT) : impl(data, deltaT) {}

    template<class ScalarType>
    FFT<ScalarType, 1>& FFT<ScalarType, 1>::operator=(FFT<ScalarType, 1> fft) noexcept {
        swap(fft);
        return *this;
    }

    template<class ScalarType>
    typename FFT<ScalarType, 1>::ComplexType FFT<ScalarType, 1>::getFreqIntense(const RealType& freq) const noexcept {
        const double float_index = double(getDeltaT() * freq * RealType(getSize()));
        double round_helper;
        if constexpr (isComplex) {
            round_helper = freq.isPositive() ? 0.5 : -0.5;
        }
        else {
            assert(!freq.isNegative());
            round_helper = 0.5;
        }
        const ssize_t index = static_cast<ssize_t>(float_index + round_helper);
        return getComponent(index);
    }

    template<class ScalarType, size_t Dim>
    class FFT {
        using Impl = Internal::FFTImpl<ScalarType, Dim>;
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
        static constexpr bool isComplex = ScalarType::isComplex;
    private:
        Impl impl;
    public:
        FFT(Utils::Array<size_t, Dim> size, Utils::Array<RealType, Dim> deltaTs);
        FFT(const Vector<ScalarType>& data, Utils::Array<size_t, Dim> size, Utils::Array<RealType, Dim> deltaTs);
        FFT(const FFT&) = default;
        FFT(FFT&&) noexcept = default;
        ~FFT() = default;
        /* Operators */
        FFT& operator=(FFT fft) noexcept;
        ComplexType operator()(size_t i) { return impl(i); }
        /* Operations */
        void transform(const Vector<ScalarType>& data) { impl.transform(data); }
        void invTransform(const Vector<ComplexType>& data) { impl.invTransform(data); }
        /* Getters */
        [[nodiscard]] size_t getSize(size_t dim) const noexcept { return impl.getSize(dim); }
        [[nodiscard]] const RealType& getDeltaT(size_t dim) const noexcept { return impl.getDeltaT(dim); }
        [[nodiscard]] ComplexType getComponent(Utils::Array<ssize_t, Dim> indexes) const { return impl.getComponent(indexes); }
        [[nodiscard]] Vector<ComplexType> getComponents() const { return impl.getComponents(); }
        [[nodiscard]] RealType getDeltaFreq(size_t dim) const noexcept { return reciprocal(impl.getDeltaT(dim) * impl.getSize(dim)); }
        [[nodiscard]] ComplexType getFreqIntense(Utils::Array<RealType, Dim> freq) const noexcept;
        /* Helpers */
        void swap(FFT& fft) noexcept { impl.swap(fft.impl); }
    };

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(Utils::Array<size_t, Dim> size, Utils::Array<RealType, Dim> deltaTs)
            : impl(size, deltaTs) {}

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>::FFT(const Vector<ScalarType>& data, Utils::Array<size_t, Dim> size, Utils::Array<RealType, Dim> deltaTs)
            : impl(data, size, deltaTs) {}

    template<class ScalarType, size_t Dim>
    FFT<ScalarType, Dim>& FFT<ScalarType, Dim>::operator=(FFT<ScalarType, Dim> fft) noexcept {
        swap(fft);
        return *this;
    }

    template<class ScalarType, size_t Dim>
    typename FFT<ScalarType, Dim>::ComplexType FFT<ScalarType, Dim>::getFreqIntense(Utils::Array<RealType, Dim> freq) const noexcept {
        Utils::Array<ssize_t, Dim> indexes{};
        for (size_t i = 0; i < Dim; ++i) {
            RealType f = freq[i];
            double round_helper;
            if constexpr (isComplex) {
                round_helper = f.isPositive() ? 0.5 : -0.5;
            }
            else {
                if (i == freq.getLength() - 1)
                    f.toAbs();
                round_helper = 0.5;
            }
            const double float_index = double(getDeltaT(i) * f * RealType(getSize(i)));
            indexes[i] = static_cast<ssize_t>(float_index + round_helper);
        }
        return getComponent(indexes);
    }
}

#include "FFTImpl.h"
