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

#include "Physica/Core/MultiPrecision/Differentiable.h"

namespace Physica::Core {
    template<class PlainScalar, DiffMode Mode>
    class FFT<Differentiable<PlainScalar, Mode, 1>, 1>
            : public FFTRSpace<FFT<Differentiable<PlainScalar, Mode, 1>, 1>, 1>
            , public FFTKSpace<FFT<Differentiable<PlainScalar, Mode, 1>, 1>, 1> {
        using ScalarType = Differentiable<PlainScalar, Mode, 1>;
        using This = FFT<ScalarType, 1>;
        using Traits = Internal::Traits<This>;
        using RealType = typename Traits::RealType;
        using ComplexType = typename Traits::ComplexType;
        using RSpaceType = FFTRSpace<This, 1>;
        using KSpaceType = FFTKSpace<This, 1>;
    private:
        using PlainRealType = typename PlainScalar::RealType;
        using PlainComplexType = typename PlainScalar::ComplexType;
        constexpr static bool isComplex = Traits::isComplex;
        using FFTType = FFT<typename std::conditional<isComplex, PlainComplexType, PlainRealType>::type, 1>;

        FFTType fft_impl;
        Utils::Array<ComplexType> buffer;
    public:
        FFT();
        FFT(size_t rSpaceSize, PlainRealType rSpaceDelta, PlanFlag planFlag);
        FFT(const Vector<ScalarType>& data, PlainRealType rSpaceDelta, PlanFlag planFlag);
        FFT(const FFT& fft) = default;
        FFT(FFT&& fft) noexcept = default;
        ~FFT() = default;
        /* Operators */
        FFT& operator=(FFT obj) noexcept;
        /* Operations */
        using RSpaceType::transform;
        using KSpaceType::invTransform;
        void transform() { transform(*this, *this); }
        void invTransform() { invTransform(*this, *this); }
        void swap(FFT& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getDim() { return 1; }
        [[nodiscard]] size_t getRSpaceSize() const noexcept { return fft_impl.getRSpaceSize(); }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return getKSpace().getSize(); }
        [[nodiscard]] RealType getRSpaceDelta() const noexcept { return fft_impl.getRSpaceDelta(); }
        [[nodiscard]] RealType getKSpaceDelta() const noexcept { return getKSpace().getDelta(); }
        [[nodiscard]] RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] inline static This makeEmptyFFT(size_t rSpaceSize, PlainRealType rSpaceDelta);
        template<class IndexType>
        [[nodiscard]] inline static IndexType rSizeToKSize(IndexType rSize) noexcept { return FFTType::rSizeToKSize(rSize); }
        inline static void transform(const This& planProvider, This& bufferProvider);
        inline static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(size_t rSpaceSize, PlainRealType rSpaceDelta);
        /* Getters */
        [[nodiscard]] RealType* asRealBuffer() { return reinterpret_cast<RealType*>(buffer.data()); }
        [[nodiscard]] const RealType* asRealBuffer() const { return reinterpret_cast<const RealType*>(buffer.data()); }
        [[nodiscard]] ComplexType* asComplexBuffer() { return buffer.data(); }
        [[nodiscard]] const ComplexType* asComplexBuffer() const { return buffer.data(); }
        /* Frients */
        friend class FFTRSpace<This, 1>;
        friend class FFTKSpace<This, 1>;
    };

    template<class PlainScalar, DiffMode Mode>
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>::FFT() : fft_impl(), buffer() {}

    template<class PlainScalar, DiffMode Mode>
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>::FFT(size_t rSpaceSize, PlainRealType rSpaceDelta, PlanFlag planFlag)
            : fft_impl(rSpaceSize, rSpaceDelta, planFlag) {
        buffer.resize(getKSpaceSize());
    }

    template<class PlainScalar, DiffMode Mode>
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>::FFT(const Vector<ScalarType>& data, PlainRealType rSpaceDelta, PlanFlag planFlag)
            : FFT(data.getLength(), rSpaceDelta, planFlag) {
        RSpaceType::transform(data);
    }

    template<class PlainScalar, DiffMode Mode>
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>::FFT(size_t rSpaceSize, PlainRealType rSpaceDelta)
            : fft_impl(FFTType::makeEmptyFFT(rSpaceSize, rSpaceDelta)) {
        buffer.resize(getKSpaceSize());
    }

    template<class PlainScalar, DiffMode Mode>
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>&
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>::operator=(FFT<Differentiable<PlainScalar, Mode, 1>, 1> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class PlainScalar, DiffMode Mode>
    void FFT<Differentiable<PlainScalar, Mode, 1>, 1>::swap(FFT& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        fft_impl.swap(obj.fft_impl);
        buffer.swap(obj.buffer);
    }

    template<class PlainScalar, DiffMode Mode>
    inline FFT<Differentiable<PlainScalar, Mode, 1>, 1>
    FFT<Differentiable<PlainScalar, Mode, 1>, 1>::makeEmptyFFT(size_t rSpaceSize, PlainRealType rSpaceDelta) {
        return This(rSpaceSize, rSpaceDelta);
    }

    template<class PlainScalar, DiffMode Mode>
    inline void FFT<Differentiable<PlainScalar, Mode, 1>, 1>::transform(const This& planProvider, This& bufferProvider) {
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        const size_t rSpaceSize = planProvider.getRSpaceSize();
        const size_t kSpaceSize = planProvider.getKSpaceSize();
        auto& rSpaceImpl = bufferProvider.fft_impl.getRSpace();
        auto& rSpace = bufferProvider.getRSpace();
        for (size_t i = 0; i < rSpaceSize; ++i)
            rSpaceImpl[i] = rSpace[i].getValue();
        FFTType::transform(planProvider.fft_impl, bufferProvider.fft_impl);

        auto& kSpaceImpl = bufferProvider.fft_impl.getKSpace();
        auto& kSpace = bufferProvider.getKSpace();
        if constexpr (isComplex) {
            for (size_t i = 0; i < rSpaceSize; ++i) {
                ComplexType& buffer_temp = rSpace[i];
                PlainComplexType& plan_temp = rSpaceImpl[i];
                buffer_temp.getValue() = plan_temp;
                plan_temp = buffer_temp.getGrad();
            }
        }
        else {
            size_t i = 0;
            for (; i < kSpaceSize - 1; ++i) {
                const PlainComplexType copy = kSpaceImpl[i];
                rSpaceImpl[2 * i] = rSpace[2 * i].getGrad();
                rSpaceImpl[2 * i + 1] = rSpace[2 * i + 1].getGrad();
                kSpace[i].getValue() = copy;
            }
            kSpace[i].getValue() = kSpaceImpl[i];
        }
        FFTType::transform(planProvider.fft_impl, bufferProvider.fft_impl);

        for (size_t i = 0; i < kSpaceSize; ++i)
            kSpace[i].getGrad() = kSpaceImpl[i];
    }

    template<class PlainScalar, DiffMode Mode>
    inline void FFT<Differentiable<PlainScalar, Mode, 1>, 1>::invTransform(const This& planProvider, This& bufferProvider) {
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        const size_t rSpaceSize = planProvider.getRSpaceSize();
        const size_t kSpaceSize = planProvider.getKSpaceSize();
        auto& kSpaceImpl = bufferProvider.fft_impl.getKSpace();
        auto& kSpace = bufferProvider.getKSpace();
        for (size_t i = 0; i < kSpaceSize; ++i)
            kSpaceImpl[i] = kSpace[i].getValue();
        FFTType::invTransform(planProvider.fft_impl, bufferProvider.fft_impl);

        auto& rSpaceImpl = bufferProvider.fft_impl.getRSpace();
        auto& rSpace = bufferProvider.getRSpace();
        if constexpr (isComplex) {
            for (size_t i = 0; i < rSpaceSize; ++i) {
                ComplexType& buffer_temp = rSpace[i];
                PlainComplexType& plan_temp = rSpaceImpl[i];
                buffer_temp.getValue() = plan_temp;
                plan_temp = buffer_temp.getGrad();
            }
        }
        else {
            size_t i = 0;
            for (; i < kSpaceSize - 1; ++i) {
                const PlainComplexType copy = kSpace[i].getGrad();
                rSpace[2 * i].getValue() = rSpaceImpl[2 * i];
                rSpace[2 * i + 1].getValue() = rSpaceImpl[2 * i + 1];
                kSpaceImpl[i] = copy;
            }
            kSpaceImpl[i] = kSpace[i].getGrad();
        }
        FFTType::invTransform(planProvider.fft_impl, bufferProvider.fft_impl);

        for (size_t i = 0; i < rSpaceSize; ++i)
            rSpace[i].getGrad() = rSpaceImpl[i];
    }
}
