/*
 * Copyright 2020-2024 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/ContinuousMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/LValueTensor.h"
#include "FFTImpl/FFTRSpace.h"
#include "FFTImpl/FFTKSpace.h"

namespace Physica {
    template<Scalar, size_t Dim = 1> class FFT;

    enum class PlanFlag {
        Measure = FFTW_MEASURE,
        Estimate = FFTW_ESTIMATE,
        Patient = FFTW_PATIENT,
        Exhaustive = FFTW_EXHAUSTIVE
    };
    /**
     * A FFT transform a tensor in r-space to a tensor in k-space.
     */
    template<Scalar T>
    class FFT<T, 1>
            : public FFTRSpace<FFT<T, 1>, 1>
            , public FFTKSpace<FFT<T, 1>, 1> {
        static_assert(!T::isDiffable, "[Error]: Header of differentiable fft should be included");
        using This = FFT<T, 1>;
        using MachineType = Traits<This>::MachineType;
        using RealType = Traits<This>::RealType;
        using ComplexType = Traits<This>::ComplexType;
        using PlanType = Traits<This>::PlanType;
        using ComplexTypeFFTW = Traits<This>::ComplexTypeFFTW;
        constexpr static bool isComplex = Traits<This>::isComplex;
        constexpr static bool isSinglePrec = Traits<This>::isSinglePrec;
    public:
        using RSpaceType = FFTRSpace<This, 1>;
        using KSpaceType = FFTKSpace<This, 1>;

        using ScalarType = RSpaceType::ScalarType;
    private:
        PlanType forward_plan;
        PlanType backward_plan;
        ComplexTypeFFTW* buffer;
        int rSpaceSize;
        PlanFlag planFlag;
    public:
        FFT();
        FFT(size_t rSpaceSize_, PlanFlag planFlag);
        FFT(const VectorND<ScalarType>& data, PlanFlag planFlag);
        FFT(const FFT& fft);
        FFT(FFT&& fft) noexcept;
        ~FFT();
        /* Operators */
        FFT& operator=(FFT obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using RSpaceType::transform;
        using KSpaceType::invTransform;
        void transform() { transform(*this, *this); }
        void rawInvTransform() { rawInvTransform(*this, *this); }
        void invTransform() { invTransform(*this, *this); }
        void swap(FFT& __restrict fft) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getDim() { return 1; }
        [[nodiscard]] __host__ __device__ ComplexTypeFFTW* getBuffer() { return buffer; }
        [[nodiscard]] __host__ __device__ const ComplexTypeFFTW* getBuffer() const { return buffer; }
        [[nodiscard]] __host__ __device__ size_t getRSpaceSize() const noexcept { return rSpaceSize; }
        [[nodiscard]] __host__ __device__ size_t getKSpaceSize() const noexcept { return getKSpace().getSize(); }
        [[nodiscard]] inline RealType getRSpaceDelta(RealType kSpaceDelta) const noexcept;
        [[nodiscard]] inline RealType getKSpaceDelta(RealType rSpaceDelta) const noexcept;
        [[nodiscard]] __host__ __device__ RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] __host__ __device__ KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] __host__ __device__ const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] __host__ __device__ const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] inline static FFT<T, 1> makeEmptyFFT(size_t rSpaceSize);
        template<class IndexType>
        [[nodiscard]] __host__ __device__ constexpr inline static IndexType rSizeToKSize(IndexType rSize) noexcept;
        static void transform(const This& planProvider, This& bufferProvider);
        static void rawInvTransform(const This& planProvider, This& bufferProvider);
        static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(size_t rSpaceSize_);
        /* Operations */
        void initializePlan();
        /* Getters */
        [[nodiscard]] __host__ __device__ RealType* asRealBuffer() { return reinterpret_cast<RealType*>(buffer); }
        [[nodiscard]] __host__ __device__ ComplexType* asComplexBuffer() { return reinterpret_cast<ComplexType*>(buffer); }
        [[nodiscard]] __host__ __device__ const RealType* asRealBuffer() const { return reinterpret_cast<const RealType*>(buffer); }
        [[nodiscard]] __host__ __device__ const ComplexType* asComplexBuffer() const { return reinterpret_cast<const ComplexType*>(buffer); }
        /* Friends */
        friend class FFTRSpace<This, 1>;
        friend class FFTKSpace<This, 1>;
    };

    template<Scalar T, size_t Dim>
    class FFT : public FFTRSpace<FFT<T, Dim>, Dim>
              , public FFTKSpace<FFT<T, Dim>, Dim> {
        static_assert(!T::isDiffable, "[Error]: Not implemented");
        using This = FFT<T, Dim>;
        using MachineType = Traits<This>::MachineType;
        using RealType = Traits<This>::RealType;
        using ComplexType = Traits<This>::ComplexType;
        using PlanType = Traits<This>::PlanType;
        using ComplexTypeFFTW = Traits<This>::ComplexTypeFFTW;
        using IndexArray = Array<size_t, Dim>;
        constexpr static bool isComplex = Traits<This>::isComplex;
        constexpr static bool isSinglePrec = Traits<This>::isSinglePrec;
    public:
        using RSpaceType = FFTRSpace<This, Dim>;
        using KSpaceType = FFTKSpace<This, Dim>;

        using ScalarType = RSpaceType::ScalarType;
    private:
        PlanType forward_plan;
        PlanType backward_plan;
        ComplexTypeFFTW* buffer;
        Array<int, Dim> rSpaceSize;
        Array<int, Dim> kSpaceSize;
        PlanFlag planFlag;
    public:
        FFT();
        FFT(const Array<size_t, Dim>& rSpaceSize_, PlanFlag planFlag_);
        FFT(const FFT&);
        FFT(FFT&&) noexcept;
        ~FFT();
        /* Operators */
        FFT& operator=(FFT obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using RSpaceType::transform;
        using KSpaceType::invTransform;
        void transform() { transform(*this, *this); }
        void rawInvTransform() { rawInvTransform(*this, *this); }
        void invTransform() { invTransform(*this, *this); }
        void swap(FFT& __restrict fft) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDim() const noexcept { return Dim == Dynamic ? rSpaceSize.getLength() : Dim; }
        [[nodiscard]] __host__ __device__ ComplexTypeFFTW* getBuffer() { return buffer; }
        [[nodiscard]] __host__ __device__ const ComplexTypeFFTW* getBuffer() const { return buffer; }
        [[nodiscard]] IndexArray getRSpaceSize() const noexcept;
        [[nodiscard]] IndexArray getKSpaceSize() const noexcept;
        [[nodiscard]] inline RealType getRSpaceDelta(RealType kSpaceDelta, unsigned int dim) const noexcept;
        [[nodiscard]] inline RealType getKSpaceDelta(RealType rSpaceDelta, unsigned int dim) const noexcept;
        [[nodiscard]] __host__ __device__ RSpaceType& getRSpace() { return *this; }
        [[nodiscard]] __host__ __device__ KSpaceType& getKSpace() { return *this; }
        [[nodiscard]] __host__ __device__ const RSpaceType& getRSpace() const { return *this; }
        [[nodiscard]] __host__ __device__ const KSpaceType& getKSpace() const { return *this; }
        /* Static members */
        [[nodiscard]] inline static FFT<T, Dim> makeEmptyFFT(const Array<size_t, Dim>& rSpaceSize);
        template<class IndexType>
        [[nodiscard]] static Array<IndexType, Dim> rSizeToKSize(const Array<IndexType, Dim>& rSize);
        static void transform(const This& planProvider, This& bufferProvider);
        static void rawInvTransform(const This& planProvider, This& bufferProvider);
        static void invTransform(const This& planProvider, This& bufferProvider);
    private:
        FFT(const Array<size_t, Dim>& rSpaceSize_);
        /* Operations */
        void initializePlan();
        PlanType makeForwardPlan();
        PlanType makeBackwardPlan();
        size_t sumRSpaceSize(size_t from_dim) const;
        size_t sumKSpaceSize(size_t from_dim) const;
        void normalizeIndexes(Array<ssize_t, Dim>& indexes) const;
        [[nodiscard]] size_t componentsSizeFrom(size_t dim) const;
        [[nodiscard]] Array<ssize_t, Dim> linearIndexToDim(size_t index) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ RealType* asRealBuffer() { return reinterpret_cast<RealType*>(buffer); }
        [[nodiscard]] __host__ __device__ ComplexType* asComplexBuffer() { return reinterpret_cast<ComplexType*>(buffer); }
        [[nodiscard]] __host__ __device__ const RealType* asRealBuffer() const { return reinterpret_cast<const RealType*>(buffer); }
        [[nodiscard]] __host__ __device__ const ComplexType* asComplexBuffer() const { return reinterpret_cast<const ComplexType*>(buffer); }
        /* Static members */
        static bool checkSize(const Array<size_t, Dim>& rSpaceSize);
        /* Friends */
        friend class FFTRSpace<This, Dim>;
        friend class FFTKSpace<This, Dim>;
    };
}

namespace Physica {
    template<Scalar T, size_t dim>
    class Traits<FFT<T, dim>> {
    public :
        using ScalarType = T;
        using MachineType = ScalarType::MachineType;
        using RealType = ScalarType::RealType;
        using ComplexType = ScalarType::ComplexType;
        constexpr static size_t Dim = dim;

        constexpr static bool isComplex = T::isComplex;
        constexpr static bool isDiffable = T::isDiffable;
        constexpr static bool isSinglePrec = std::is_same<MachineType, float>::value;
        using PlanType = std::conditional<isSinglePrec, fftwf_plan, fftw_plan>::type;
        using ComplexTypeFFTW = std::conditional<isSinglePrec, fftwf_complex, fftw_complex>::type;

        static_assert(isDiffable || sizeof(RealType) == sizeof(MachineType), "[Error]: Invalid ScalarType");
        static_assert(isDiffable || sizeof(ComplexType) == sizeof(ComplexTypeFFTW), "[Error]: Invalid ScalarType");
        static_assert(Dim <= 3U, "[Error]: Dimension higher than 3 should be declared as dynamic");
        static_assert(Dim != 0, "[Error]: Not implemented");
    };
}

#include "FFTImpl/FFTImpl.h"
