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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/LValueGrid.h"
#include "FFTImpl/FFTRSpace.h"
#include "FFTImpl/FFTKSpace.h"

namespace Physica::Core {
    template<class ScalarType, size_t Dim = 1> class FFT;

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
            constexpr static bool isSinglePrec = std::is_same<TrivialType, float>::value;
            using PlanType = typename std::conditional<isSinglePrec, fftwf_plan, fftw_plan>::type;
            using ComplexTypeFFTW = typename std::conditional<isSinglePrec, fftwf_complex, fftw_complex>::type;

            static_assert(sizeof(RealType) == sizeof(TrivialType), "[Error]: Invalid ScalarType");
            static_assert(sizeof(ComplexType) == sizeof(ComplexTypeFFTW), "[Error]: Invalid ScalarType");

            enum PlanFlag {
                Measure = FFTW_MEASURE,
                Estimate = FFTW_ESTIMATE,
                Patient = FFTW_PATIENT,
                Exhaustive = FFTW_EXHAUSTIVE
            };
        };
    }
    /**
     * A FFT transform a tensor in r-space to a tensor in k-space.
     */
    template<class ScalarType>
    class FFT<ScalarType, 1>
            : public Internal::Traits<FFT<ScalarType, 1>>
            , public FFTRSpace<FFT<ScalarType, 1>, 1>
            , public FFTKSpace<FFT<ScalarType, 1>, 1> {
        using This = FFT<ScalarType, 1>;
        using Traits = Internal::Traits<This>;
        using typename Traits::TrivialType;
        using typename Traits::RealType;
        using typename Traits::ComplexType;
        using typename Traits::PlanType;
        using typename Traits::ComplexTypeFFTW;
        using Traits::isComplex;
        using Traits::isSinglePrec;
        using RSpaceType = FFTRSpace<This, 1>;
        using KSpaceType = FFTKSpace<This, 1>;
    public:
        using typename Traits::PlanFlag;
    private:
        PlanType forward_plan;
        PlanType backward_plan;
        union {
            RealType* real_buffer;
            ComplexType* complex_buffer;
            ComplexTypeFFTW* buffer;
        };
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
        inline void transform();
        inline void invTransform();
        void swap(FFT& fft) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getDim() { return 1; }
        [[nodiscard]] size_t getRSpaceSize() const noexcept { return getRSpace().getSize(); }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return getKSpace().getSize(); }
        [[nodiscard]] const RealType& getRSpaceDelta() const noexcept { return rSpaceDelta; }
        [[nodiscard]] RealType getKSpaceDelta() const noexcept { return reciprocal(getRSpaceDelta() * getRSpaceSize()); }
        [[nodiscard]] RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] inline static FFT<ScalarType, 1> makeEmptyFFT(size_t rSpaceSize, RealType rSpaceDelta);
        [[nodiscard]] inline static int rSizeToKSize(int rSize) noexcept;
        static void transform(const This& planProvider, This& bufferProvider);
        static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(size_t rSpaceSize_, RealType rSpaceDelta_);
        /* Operations */
        void initializePlan();
        /* Friends */
        friend class FFTRSpace<This, 1>;
        friend class FFTKSpace<This, 1>;
    };

    template<class ScalarType, size_t Dim>
    class FFT
            : public Internal::Traits<FFT<ScalarType, Dim>>
            , public FFTRSpace<FFT<ScalarType, Dim>, Dim>
            , public FFTKSpace<FFT<ScalarType, Dim>, Dim> {
        static_assert(Dim <= 3U, "[Error]: Dimension higher than 3 should be declared as dynamic");
        static_assert(Dim != 0, "[Error]: Not implemented");
        using This = FFT<ScalarType, Dim>;
        using Traits = Internal::Traits<This>;
        using typename Traits::TrivialType;
        using typename Traits::RealType;
        using typename Traits::ComplexType;
        using typename Traits::PlanType;
        using typename Traits::ComplexTypeFFTW;
        using Traits::isComplex;
        using Traits::isSinglePrec;
        using RSpaceType = FFTRSpace<This, Dim>;
        using KSpaceType = FFTKSpace<This, Dim>;
    public:
        using typename Traits::PlanFlag;
    private:
        PlanType forward_plan;
        PlanType backward_plan;
        union {
            RealType* real_buffer;
            ComplexType* complex_buffer;
            ComplexTypeFFTW* buffer;
        };
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
        void transform();
        void invTransform();
        void swap(FFT& fft) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return Dim == Dynamic ? rSpaceSize.getLength() : Dim; }
        [[nodiscard]] const Utils::Array<int, Dim>& getRSpaceSize() const noexcept { return rSpaceSize; }
        [[nodiscard]] const Utils::Array<int, Dim>& getKSpaceSize() const noexcept { return kSpaceSize; }
        [[nodiscard]] RealType getRSpaceDelta(size_t dim) const noexcept { return rSpaceDelta[dim]; }
        [[nodiscard]] RealType getKSpaceDelta(size_t dim) const noexcept { return reciprocal(getRSpaceDelta(dim) * getRSpaceSize()[dim]); }
        [[nodiscard]] RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] static Utils::Array<int, Dim> rSizeToKSize(const Utils::Array<int, Dim>& rSize);
        static void transform(const This& planProvider, This& bufferProvider);
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
        /* Static members */
        static bool checkSize(const Utils::Array<size_t, Dim>& rSpaceSize);
        /* Friends */
        friend class FFTRSpace<This, Dim>;
        friend class FFTKSpace<This, Dim>;
    };
}

#include "FFTImpl/FFTImpl.h"
