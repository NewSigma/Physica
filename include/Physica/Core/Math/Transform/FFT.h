/*
 * Copyright 2020-2023 WeiBo He.
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

#include <fftw3.h>
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/ContinuousMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/GridImpl/LValueGrid.h"
#include "FFTImpl/FFTRSpace.h"
#include "FFTImpl/FFTKSpace.h"

namespace Physica::Core {
    template<class ScalarType, size_t Dim = 1> class FFT;

    enum class PlanFlag {
        Measure = FFTW_MEASURE,
        Estimate = FFTW_ESTIMATE,
        Patient = FFTW_PATIENT,
        Exhaustive = FFTW_EXHAUSTIVE
    };

    namespace Internal {
        template<class T, size_t dim>
        class Traits<FFT<T, dim>> {
        public :
            using ScalarType = T;
            using TrivialType = typename ScalarType::TrivialType;
            using RealType = typename ScalarType::RealType;
            using ComplexType = typename ScalarType::ComplexType;
            constexpr static size_t Dim = dim;

            constexpr static bool isComplex = T::isComplex;
            constexpr static bool isDifferentiable = T::isDifferentiable;
            constexpr static bool isSinglePrec = std::is_same<TrivialType, float>::value;
            using PlanType = typename std::conditional<isSinglePrec, fftwf_plan, fftw_plan>::type;
            using ComplexTypeFFTW = typename std::conditional<isSinglePrec, fftwf_complex, fftw_complex>::type;

            static_assert(isDifferentiable || sizeof(RealType) == sizeof(TrivialType), "[Error]: Invalid ScalarType");
            static_assert(isDifferentiable || sizeof(ComplexType) == sizeof(ComplexTypeFFTW), "[Error]: Invalid ScalarType");
        };
    }
    /**
     * A FFT transform a tensor in r-space to a tensor in k-space.
     */
    template<class ScalarType>
    class FFT<ScalarType, 1>
            : public FFTRSpace<FFT<ScalarType, 1>, 1>
            , public FFTKSpace<FFT<ScalarType, 1>, 1> {
        using This = FFT<ScalarType, 1>;
        using Traits = Internal::Traits<This>;
        using TrivialType = typename Traits::TrivialType;
        using RealType = typename Traits::RealType;
        using ComplexType = typename Traits::ComplexType;
        using PlanType = typename Traits::PlanType;
        using ComplexTypeFFTW = typename Traits::ComplexTypeFFTW;
        constexpr static bool isComplex = Traits::isComplex;
        constexpr static bool isSinglePrec = Traits::isSinglePrec;
        static_assert(!Traits::isDifferentiable, "[Error]: Header of differentiable fft should be included");
    public:
        using RSpaceType = FFTRSpace<This, 1>;
        using KSpaceType = FFTKSpace<This, 1>;
    private:
        PlanType forward_plan;
        PlanType backward_plan;
        ComplexTypeFFTW* buffer;
        int rSpaceSize;
        RealType rSpaceDelta;
        PlanFlag planFlag;
    public:
        FFT();
        FFT(size_t rSpaceSize_, RealType rSpaceDelta_, PlanFlag planFlag);
        FFT(const Vector<ScalarType>& data, RealType rSpaceDelta_, PlanFlag planFlag);
        FFT(const FFT& fft);
        FFT(FFT&& fft) noexcept;
        ~FFT();
        /* Operators */
        FFT& operator=(FFT fft) noexcept;
        /* Operations */
        using RSpaceType::transform;
        using KSpaceType::invTransform;
        void transform() { transform(*this, *this); }
        void rawInvTransform() { rawInvTransform(*this, *this); }
        void invTransform() { invTransform(*this, *this); }
        void swap(FFT& fft) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getDim() { return 1; }
        [[nodiscard]] ComplexTypeFFTW* getBuffer() { return buffer; }
        [[nodiscard]] const ComplexTypeFFTW* getBuffer() const { return buffer; }
        [[nodiscard]] size_t getRSpaceSize() const noexcept { return rSpaceSize; }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return getKSpace().getSize(); }
        [[nodiscard]] RealType getRSpaceDelta() const noexcept { return rSpaceDelta; }
        [[nodiscard]] RealType getKSpaceDelta() const noexcept { return getKSpace().getDelta(); }
        [[nodiscard]] RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] inline static FFT<ScalarType, 1> makeEmptyFFT(size_t rSpaceSize, RealType rSpaceDelta);
        template<class IndexType>
        [[nodiscard]] inline static IndexType rSizeToKSize(IndexType rSize) noexcept;
        static void transform(const This& planProvider, This& bufferProvider);
        static void rawInvTransform(const This& planProvider, This& bufferProvider);
        static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(size_t rSpaceSize_, RealType rSpaceDelta_);
        /* Operations */
        void initializePlan();
        /* Getters */
        [[nodiscard]] RealType* asRealBuffer() { return reinterpret_cast<RealType*>(buffer); }
        [[nodiscard]] ComplexType* asComplexBuffer() { return reinterpret_cast<ComplexType*>(buffer); }
        [[nodiscard]] const RealType* asRealBuffer() const { return reinterpret_cast<const RealType*>(buffer); }
        [[nodiscard]] const ComplexType* asComplexBuffer() const { return reinterpret_cast<const ComplexType*>(buffer); }
        /* Friends */
        friend class FFTRSpace<This, 1>;
        friend class FFTKSpace<This, 1>;
    };

    template<class ScalarType, size_t Dim>
    class FFT : public FFTRSpace<FFT<ScalarType, Dim>, Dim>
              , public FFTKSpace<FFT<ScalarType, Dim>, Dim> {
        static_assert(Dim <= 3U, "[Error]: Dimension higher than 3 should be declared as dynamic");
        static_assert(Dim != 0, "[Error]: Not implemented");
        static_assert(!ScalarType::isDifferentiable, "[Error]: Not implemented");
        using This = FFT<ScalarType, Dim>;
        using Traits = Internal::Traits<This>;
        using TrivialType = typename Traits::TrivialType;
        using RealType = typename Traits::RealType;
        using ComplexType = typename Traits::ComplexType;
        using PlanType = typename Traits::PlanType;
        using ComplexTypeFFTW = typename Traits::ComplexTypeFFTW;
        using IndexArray = Utils::Array<size_t, Dim>;
        constexpr static bool isComplex = Traits::isComplex;
        constexpr static bool isSinglePrec = Traits::isSinglePrec;
    public:
        using RSpaceType = FFTRSpace<This, Dim>;
        using KSpaceType = FFTKSpace<This, Dim>;
    private:
        PlanType forward_plan;
        PlanType backward_plan;
        ComplexTypeFFTW* buffer;
        Utils::Array<int, Dim> rSpaceSize;
        Utils::Array<int, Dim> kSpaceSize;
        Utils::Array<RealType, Dim> rSpaceDelta;
        PlanFlag planFlag;
    public:
        FFT();
        FFT(const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_, PlanFlag planFlag_);
        FFT(const FFT&);
        FFT(FFT&&) noexcept;
        ~FFT();
        /* Operators */
        FFT& operator=(FFT fft) noexcept;
        /* Operations */
        using RSpaceType::transform;
        using KSpaceType::invTransform;
        void transform() { transform(*this, *this); }
        void rawInvTransform() { rawInvTransform(*this, *this); }
        void invTransform() { invTransform(*this, *this); }
        void swap(FFT& fft) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return Dim == Dynamic ? rSpaceSize.getLength() : Dim; }
        [[nodiscard]] ComplexTypeFFTW* getBuffer() { return buffer; }
        [[nodiscard]] const ComplexTypeFFTW* getBuffer() const { return buffer; }
        [[nodiscard]] IndexArray getRSpaceSize() const noexcept;
        [[nodiscard]] IndexArray getKSpaceSize() const noexcept;
        [[nodiscard]] RealType getRSpaceDelta(size_t dim) const noexcept { return rSpaceDelta[dim]; }
        [[nodiscard]] RealType getKSpaceDelta(size_t dim) const noexcept { return reciprocal(getRSpaceDelta(dim) * getRSpaceSize()[dim]); }
        [[nodiscard]] RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        template<class IndexType>
        [[nodiscard]] static Utils::Array<IndexType, Dim> rSizeToKSize(const Utils::Array<IndexType, Dim>& rSize);
        static void transform(const This& planProvider, This& bufferProvider);
        static void rawInvTransform(const This& planProvider, This& bufferProvider);
        static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_);
        /* Operations */
        void initializePlan();
        fftw_plan makeForwardPlan();
        fftw_plan makeBackwardPlan();
        size_t sumRSpaceSize(size_t from_dim) const;
        size_t sumKSpaceSize(size_t from_dim) const;
        void normalizeIndexes(Utils::Array<ssize_t, Dim>& indexes) const;
        [[nodiscard]] RealType mulDeltas() const;
        [[nodiscard]] size_t componentsSizeFrom(size_t dim) const;
        [[nodiscard]] Utils::Array<ssize_t, Dim> linearIndexToDim(size_t index) const;
        /* Getters */
        [[nodiscard]] RealType* asRealBuffer() { return reinterpret_cast<RealType*>(buffer); }
        [[nodiscard]] ComplexType* asComplexBuffer() { return reinterpret_cast<ComplexType*>(buffer); }
        [[nodiscard]] const RealType* asRealBuffer() const { return reinterpret_cast<const RealType*>(buffer); }
        [[nodiscard]] const ComplexType* asComplexBuffer() const { return reinterpret_cast<const ComplexType*>(buffer); }
        /* Static members */
        static bool checkSize(const Utils::Array<size_t, Dim>& rSpaceSize);
        /* Friends */
        friend class FFTRSpace<This, Dim>;
        friend class FFTKSpace<This, Dim>;
    };
}

#include "FFTImpl/FFTImpl.h"
#include "FFTImpl/DifferentiableFFT.h"
