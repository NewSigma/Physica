/*
 * Copyright 2024 Weibo He.
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

#include <iomanip>

namespace Physica {
    template<>
    class Traits<Core::float32> {
    public:
        constexpr static Core::ScalarOption Option = Core::Float32;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = Core::float32;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Core::Complex<ScalarType>;
        using MachineType = float;
    };
}

namespace Physica::Core {
    template<>
    class Real<Float32> : public ScalarBase<Real<Float32>> {
        using This = Real<Float32>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        float f;
    public:
        Real() = default;
        __host__ __device__ Real(float f_) : f(f_) {}
        Real(const Integer& i) : Real(float(double(i))) {}
        Real(const Rational& r) : Real(float(double(r))) {}
        template<class OtherScalar>
        __host__ __device__ explicit inline Real(const ScalarBase<OtherScalar>& s);
        Real(const Real&) = default;
        Real(Real&&) noexcept = default;
        //~Real() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        using Base::operator>;
        using Base::operator<;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const { return f; }
        __host__ __device__ explicit operator double() const { return f; }
        __host__ __device__ Real operator+(const Real& s) const { return Real(f + s.f); }
        __host__ __device__ Real operator-(const Real& s) const { return Real(f - s.f); }
        __host__ __device__ Real operator*(const Real& s) const { return Real(f * s.f); }
        __host__ __device__ Real operator/(const Real& s) const { return Real(f / s.f); }
        Real operator<<(int i) const { return Real(f * std::pow(2, i)); }
        Real operator>>(int i) const { return Real(f / std::pow(2, i)); }
        __host__ __device__ Real operator-() const noexcept { return Real(-f); }
        __host__ __device__ bool operator>(const Real& s) const { return f > s.f; }
        __host__ __device__ bool operator<(const Real& s) const { return f < s.f; }
        __host__ __device__ bool operator==(const Real& s) const { return f == s.f; }
        PHYSICA_API friend std::istream& operator>>(std::istream& is, Real& scalar);
        /* Operations */
        void swap(Real& __restrict s) noexcept { std::swap(f, s.f); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float32; }
        [[nodiscard]] __host__ __device__ float toMachine() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return f == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return f > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return f < 0; }
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        [[nodiscard]] bool isInteger() const;
        /* Static Members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Real random_uniform(RandomGenerator& gen);
        template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
        [[nodiscard]] inline static Real random_uniform(Random<RandomGenerator, FixedSeed>& pool);
        template<class RandomGenerator>
        [[nodiscard]] inline static Real random_normal(RandomGenerator& gen);
        template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
        [[nodiscard]] inline static Real random_normal(Random<RandomGenerator, FixedSeed>& pool);
        template<class RandomType>
        [[nodiscard]] static Real random_normal(GaussRandomPool<This, RandomType>& pool) { return pool(); }
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Real random_any(Distribution& dist, RandomGenerator& gen);
    #ifdef PHYSICA_HDF5
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_FLOAT; }
    #endif
    };

    template<class OtherScalar>
    __host__ __device__ inline Real<Float32>::Real(const ScalarBase<OtherScalar>& s) : f(float(s.getDerived())) {}

    __host__ __device__ inline bool Real<Float32>::isFinite() const noexcept {
    #ifdef __CUDA_ARCH__
        return isfinite(f);
    #else
        return std::isfinite(f);
    #endif
    }

    template<class RandomGenerator>
    inline Real<Float32> Real<Float32>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<float> dist{};
        return Real(dist(gen));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline Real<Float32> Real<Float32>::random_uniform(Random<RandomGenerator, FixedSeed>& pool) {
        return random_uniform<RandomGenerator>(pool.getGen());
    }

    template<class RandomGenerator>
    inline Real<Float32> Real<Float32>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<float> dist{};
        return Real(dist(gen));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline Real<Float32> Real<Float32>::random_normal(Random<RandomGenerator, FixedSeed>& pool) {
        return random_normal<RandomGenerator>(pool.getGen());
    }

    template<class Distribution, class RandomGenerator>
    inline Real<Float32> Real<Float32>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Real(dist(gen));
    }

    inline std::ostream& operator<<(std::ostream& os, const Real<Float32>& s) {
        const auto lastPrec = os.precision();
        return os << std::setprecision(7) << float(s) << std::setprecision(lastPrec);
    }

    inline Real<Float32>& operator++(Real<Float32>& s) {
        s += float32(1.0F);
        return s;
    }

    inline Real<Float32>& operator--(Real<Float32>& s) {
        s -= float32(1.0F);
        return s;
    }

    inline Real<Float32> operator++(Real<Float32>& s, int) {
        float32 temp(s);
        s += float32(1.0F);
        return temp;
    }

    inline Real<Float32> operator--(Real<Float32>& s, int) {
        float32 temp(s);
        s -= float32(1.0F);
        return temp;
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Real<Physica::Core::Float>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<float>
    #else
        public numeric_limits<float>
    #endif
    {};
}
