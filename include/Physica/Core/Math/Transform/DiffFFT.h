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

#include "Physica/Core/Scalar/Diff.h"
#include "FFT.h"

namespace Physica {
    template<Scalar T, DiffMode Mode>
    class FFT<Diff<T, Mode, 1>, 1>
            : public FFTRSpace<FFT<Diff<T, Mode, 1>, 1>, 1>
            , public FFTKSpace<FFT<Diff<T, Mode, 1>, 1>, 1> {
        using ScalarType = Diff<T, Mode, 1>;
        using This = FFT<ScalarType, 1>;
        using RealType = Traits<This>::RealType;
        using RealPtrTy = RealType::PtrTy;
        using RealConstPtrTy = RealType::ConstPtrTy;
        using ComplexType = Traits<This>::ComplexType;
        using ComplexPtrTy = ComplexType::PtrTy;
        using ComplexConstPtrTy = ComplexType::ConstPtrTy;
        using RSpaceType = FFTRSpace<This, 1>;
        using KSpaceType = FFTKSpace<This, 1>;
    private:
        using PlainRealType = T::RealType;
        using PlainComplexType = T::ComplexType;
        constexpr static bool isComplex = Traits<This>::isComplex;
        using FFTType = FFT<typename std::conditional<isComplex, PlainComplexType, PlainRealType>::type, 1>;

        FFTType fft_impl;
        VectorND<ComplexType> buffer;
    public:
        FFT();
        FFT(size_t rSpaceSize, PlanFlag planFlag);
        FFT(const VectorND<ScalarType>& data, PlanFlag planFlag);
        FFT(const FFT& fft) = default;
        FFT(FFT&& fft) noexcept = default;
        ~FFT() = default;
        /* Operators */
        FFT& operator=(FFT obj) noexcept { swap(obj); return *this; }
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
        [[nodiscard]] static This makeEmptyFFT(size_t rSpaceSize);
        template<std::integral IndexType>
        [[nodiscard]] static IndexType rSizeToKSize(IndexType rSize) noexcept { return FFTType::rSizeToKSize(rSize); }
        static void transform(const This& planProvider, This& bufferProvider);
        static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(size_t rSpaceSize);
        /* Getters */
        [[nodiscard]] RealPtrTy asRealBuffer();
        [[nodiscard]] RealConstPtrTy asRealBuffer() const;
        [[nodiscard]] ComplexPtrTy asComplexBuffer();
        [[nodiscard]] ComplexConstPtrTy asComplexBuffer() const;
        /* Frients */
        friend class FFTRSpace<This, 1>;
        friend class FFTKSpace<This, 1>;
    };

    template<Scalar T, DiffMode Mode>
    FFT<Diff<T, Mode, 1>, 1>::FFT() : fft_impl(), buffer() {}

    template<Scalar T, DiffMode Mode>
    FFT<Diff<T, Mode, 1>, 1>::FFT(size_t rSpaceSize, PlanFlag planFlag)
            : fft_impl(rSpaceSize, planFlag) {
        buffer.resize(getKSpaceSize());
    }

    template<Scalar T, DiffMode Mode>
    FFT<Diff<T, Mode, 1>, 1>::FFT(const VectorND<ScalarType>& data, PlanFlag planFlag)
            : FFT(data.getLength(), planFlag) {
        RSpaceType::transform(data);
    }

    template<Scalar T, DiffMode Mode>
    FFT<Diff<T, Mode, 1>, 1>::FFT(size_t rSpaceSize)
            : fft_impl(FFTType::makeEmptyFFT(rSpaceSize)) {
        buffer.resize(getKSpaceSize());
    }

    template<Scalar T, DiffMode Mode>
    void FFT<Diff<T, Mode, 1>, 1>::swap(FFT& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        fft_impl.swap(obj.fft_impl);
        buffer.swap(obj.buffer);
    }

    template<Scalar T, DiffMode Mode>
    auto FFT<Diff<T, Mode, 1>, 1>::makeEmptyFFT(size_t rSpaceSize) -> This {
        return This(rSpaceSize);
    }

    template<Scalar T, DiffMode Mode>
    void FFT<Diff<T, Mode, 1>, 1>::transform(const This& planProvider, This& bufferProvider) {
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        const size_t rSpaceSize = planProvider.getRSpaceSize();
        const size_t kSpaceSize = planProvider.getKSpaceSize();
        auto& rSpaceImpl = bufferProvider.fft_impl.getRSpace();
        auto& rSpace = bufferProvider.getRSpace();
        for (size_t i = 0; i < rSpaceSize; ++i)
            rSpaceImpl[i] = rSpace[i].value();
        FFTType::transform(planProvider.fft_impl, bufferProvider.fft_impl);

        auto& kSpaceImpl = bufferProvider.fft_impl.getKSpace();
        auto& kSpace = bufferProvider.getKSpace();
        if constexpr (isComplex) {
            for (size_t i = 0; i < rSpaceSize; ++i) {
                ComplexType& buffer_temp = rSpace[i];
                PlainComplexType& plan_temp = rSpaceImpl[i];
                buffer_temp.value() = plan_temp;
                plan_temp = buffer_temp.grad();
            }
        }
        else {
            size_t i = 0;
            for (; i < kSpaceSize - 1; ++i) {
                const PlainComplexType copy = kSpaceImpl[i];
                rSpaceImpl[2 * i] = rSpace[2 * i].grad();
                rSpaceImpl[2 * i + 1] = rSpace[2 * i + 1].grad();
                kSpace[i].value() = copy;
            }
            kSpace[i].value() = kSpaceImpl[i];
        }
        FFTType::transform(planProvider.fft_impl, bufferProvider.fft_impl);

        for (size_t i = 0; i < kSpaceSize; ++i)
            kSpace[i].grad() = kSpaceImpl[i];
    }

    template<Scalar T, DiffMode Mode>
    void FFT<Diff<T, Mode, 1>, 1>::invTransform(const This& planProvider, This& bufferProvider) {
        assert(planProvider.getRSpaceSize() == bufferProvider.getRSpaceSize());
        assert(planProvider.getKSpaceSize() == bufferProvider.getKSpaceSize());
        const size_t rSpaceSize = planProvider.getRSpaceSize();
        const size_t kSpaceSize = planProvider.getKSpaceSize();
        auto& kSpaceImpl = bufferProvider.fft_impl.getKSpace();
        auto& kSpace = bufferProvider.getKSpace();
        for (size_t i = 0; i < kSpaceSize; ++i)
            kSpaceImpl[i] = kSpace[i].value();
        FFTType::invTransform(planProvider.fft_impl, bufferProvider.fft_impl);

        auto& rSpaceImpl = bufferProvider.fft_impl.getRSpace();
        auto& rSpace = bufferProvider.getRSpace();
        if constexpr (isComplex) {
            for (size_t i = 0; i < rSpaceSize; ++i) {
                ComplexType& buffer_temp = rSpace[i];
                PlainComplexType& plan_temp = rSpaceImpl[i];
                buffer_temp.value() = plan_temp;
                plan_temp = buffer_temp.grad();
            }
        }
        else {
            size_t i = 0;
            for (; i < kSpaceSize - 1; ++i) {
                const PlainComplexType copy = kSpace[i].grad();
                rSpace[2 * i].value() = rSpaceImpl[2 * i];
                rSpace[2 * i + 1].value() = rSpaceImpl[2 * i + 1];
                kSpaceImpl[i] = copy;
            }
            kSpaceImpl[i] = kSpace[i].grad();
        }
        FFTType::invTransform(planProvider.fft_impl, bufferProvider.fft_impl);

        for (size_t i = 0; i < rSpaceSize; ++i)
            rSpace[i].grad() = rSpaceImpl[i];
    }

    template<Scalar T, DiffMode Mode>
    auto FFT<Diff<T, Mode, 1>, 1>::asRealBuffer() -> RealPtrTy {
        return RealPtrTy(buffer.data());
    }

    template<Scalar T, DiffMode Mode>
    auto FFT<Diff<T, Mode, 1>, 1>::asRealBuffer() const -> RealConstPtrTy {
        return const_cast<This&>(*this).asRealBuffer();
    }

    template<Scalar T, DiffMode Mode>
    auto FFT<Diff<T, Mode, 1>, 1>::asComplexBuffer() -> ComplexPtrTy {
        return buffer.data();
    }

    template<Scalar T, DiffMode Mode>
    auto FFT<Diff<T, Mode, 1>, 1>::asComplexBuffer() const -> ComplexConstPtrTy {
        return const_cast<This&>(*this).asComplexBuffer();
    }
}
