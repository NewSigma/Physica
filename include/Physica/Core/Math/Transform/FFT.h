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
#include "FFTImpl/FFTRSpace.h"
#include "FFTImpl/FFTKSpace.h"

namespace Physica::Core {
    template<class ScalarType, size_t Dim = 1> class FFT;

    namespace Internal {
        template<class T, size_t dim>
        class Traits<FFT<T, dim>> {
        public :
            using ScalarType = T;
            constexpr static size_t Dim = dim;
            constexpr static bool isComplex = T::isComplex;
        };
    }
    /**
     * A FFT transform a vector in r-space to a vector in k-space.
     */
    template<class ScalarType>
    class FFT<ScalarType, 1> : public FFTRSpace<FFT<ScalarType, 1>>, public FFTKSpace<FFT<ScalarType, 1>> {
        using This = FFT<ScalarType, 1>;
        using TrivialType = typename ScalarType::TrivialType;
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;

        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static bool isSinglePrec = std::is_same<TrivialType, float>::value;
        using PlanType = typename std::conditional<isSinglePrec, fftwf_plan, fftw_plan>::type;
        using ComplexTypeFFTW = typename std::conditional<isSinglePrec, fftwf_complex, fftw_complex>::type;

        static_assert(sizeof(RealType) == sizeof(TrivialType), "[Error]: Invalid ScalarType");
        static_assert(sizeof(ComplexType) == sizeof(ComplexTypeFFTW), "[Error]: Invalid ScalarType");
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
    public:
        FFT();
        FFT(size_t rSpaceSize_, const RealType& rSpaceDelta_);
        FFT(const Vector<ScalarType>& data, const RealType& rSpaceDelta_);
        FFT(const FFT& fft);
        FFT(FFT&& fft) noexcept;
        ~FFT();
        /* Operators */
        FFT& operator=(FFT fft) noexcept;
        /* Operations */
        inline void transform();
        template<class VectorType> inline void transform(const RValueVector<VectorType>& data);
        void invTransform();
        template<class VectorType> inline void invTransform(const RValueVector<VectorType>& data);
        void swap(FFT& fft) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getDimen() { return 1; }
        [[nodiscard]] size_t getRSpaceSize() const noexcept { return getRSpace().getSize(); }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return getKSpace().getSize(); }
        [[nodiscard]] const RealType& getRSpaceDelta() const noexcept { return rSpaceDelta; }
        [[nodiscard]] RealType getKSpaceDelta() const noexcept { return reciprocal(getRSpaceDelta() * getRSpaceSize()); }
        [[nodiscard]] FFTRSpace<This>& getRSpace() { return *this; }
        [[nodiscard]] FFTKSpace<This>& getKSpace() { return *this; }
        [[nodiscard]] const FFTRSpace<This>& getRSpace() const { return *this; }
        [[nodiscard]] const FFTKSpace<This>& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] inline static int rSizeToKSize(int rSize) noexcept;
    private:
        void initializePlan();

        friend class FFTRSpace<This>;
        friend class FFTKSpace<This>;
    };

    template<class ScalarType, size_t Dim>
    class FFT {
        using RealType = typename ScalarType::RealType;
        using ComplexType = typename ScalarType::ComplexType;
        static constexpr bool isComplex = ScalarType::isComplex;
    private:
        fftw_plan forward_plan;
        fftw_plan backward_plan;
        union {
            double* real_buffer;
            fftw_complex* buffer;
        };
        Utils::Array<int, Dim> rSpaceSize;
        Utils::Array<int, Dim> kSpaceSize;
        Utils::Array<RealType, Dim> rSpaceDelta;
    public:
        FFT();
        FFT(const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_);
        FFT(const Vector<ScalarType>& data, const Utils::Array<size_t, Dim>& rSpaceSize_, Utils::Array<RealType, Dim> rSpaceDelta_);
        FFT(const FFT&);
        FFT(FFT&&) noexcept;
        ~FFT();
        /* Operators */
        FFT& operator=(FFT fft) noexcept;
        /* Operations */
        void transform(const Vector<ScalarType>& data);
        void invTransform(const Vector<ComplexType>& data);
        void swap(FFT& fft) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getDimen() { return Dim; }
        [[nodiscard]] const Utils::Array<int, Dim>& getRSpaceSize() const noexcept { return rSpaceSize; }
        [[nodiscard]] int getRSpaceSize(size_t dim) const noexcept { return rSpaceSize[dim]; }
        [[nodiscard]] const Utils::Array<int, Dim>& getKSpaceSize() const noexcept { return kSpaceSize; }
        [[nodiscard]] int getKSpaceSize(size_t dim) const noexcept { return kSpaceSize[dim]; }
        [[nodiscard]] const Utils::Array<RealType, Dim>& getRSpaceDelta() const noexcept { return rSpaceDelta; }
        [[nodiscard]] RealType getRSpaceDelta(size_t dim) const noexcept { return rSpaceDelta[dim]; }
        [[nodiscard]] RealType getKSpaceDelta(size_t dim) const noexcept { return reciprocal(getRSpaceDelta(dim) * getRSpaceSize(dim)); }
        [[nodiscard]] Vector<ScalarType> getRSpace() const;
        [[nodiscard]] Vector<ComplexType> getKSpace() const;
    private:
        static bool checkSize(const Utils::Array<size_t, Dim>& rSpaceSize);
        Utils::Array<int, Dim> makeKSpaceSize() const;
        fftw_plan forwardPlan();
        fftw_plan backwardPlan();
        void transform();
        void invTransform();
        size_t totalRSpaceSize(size_t from_dim) const;
        size_t totalKSpaceSize(size_t from_dim) const;
        void normalizeIndexes(Utils::Array<ssize_t, Dim>& indexes) const;
        [[nodiscard]] RealType mulDeltas() const;
        [[nodiscard]] size_t componentsSizeFrom(size_t dim) const;
        [[nodiscard]] Utils::Array<ssize_t, Dim> linearIndexToDim(size_t index) const;
    };
}

#include "FFTImpl/FFTImpl.h"
