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
    class Traits<Core::float64> {
    public:
        using ScalarType = Core::float64;
        using RealType = ScalarType;
        using ComplexType = Core::Complex<ScalarType>;
        using TrivialType = double;
        using PlainScalar = ScalarType;
        constexpr static Core::ScalarOption Option = Core::Float64;
        constexpr static bool isComplex = false;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;
        constexpr static unsigned int Order = 0;
    };
}

namespace Physica::Core {
    template<>
    class Scalar<Float64> : public ScalarBase<Scalar<Float64>> {
        using This = Scalar<Float64>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        double d;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(double d_) noexcept : d(d_) {}
        Scalar(const Integer& i) : Scalar(double(i)) {}
        Scalar(const Rational& r) : Scalar(double(r)) {}
        template<class OtherScalar>
        __host__ __device__ explicit inline Scalar(const ScalarBase<OtherScalar>& s);
        Scalar(const Scalar&) = default;
        Scalar(Scalar&&) noexcept = default;
        //~Scalar() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        using Base::operator>;
        using Base::operator<;
        Scalar& operator=(const Scalar& obj) = default;
        Scalar& operator=(Scalar&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const noexcept { return d; }
        __host__ __device__ explicit operator double() const noexcept { return d; }
        __host__ __device__ Scalar operator+(const Scalar& s) const noexcept { return Scalar(d + s.d); }
        __host__ __device__ Scalar operator-(const Scalar& s) const noexcept { return Scalar(d - s.d); }
        __host__ __device__ Scalar operator*(const Scalar& s) const noexcept { return Scalar(d * s.d); }
        __host__ __device__ Scalar operator/(const Scalar& s) const noexcept { return Scalar(d / s.d); }
        Scalar operator<<(int i) const { return Scalar(d * std::pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(d / std::pow(2, i)); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-d); }
        __host__ __device__ bool operator>(const Scalar& s) const noexcept { return d > s.d; }
        __host__ __device__ bool operator<(const Scalar& s) const noexcept { return d < s.d; }
        __host__ __device__ bool operator==(const Scalar& s) const noexcept { return d == s.d; }
        PHYSICA_API friend std::istream& operator>>(std::istream& is, Scalar& scalar);
        /* Operations */
        void swap(Scalar& __restrict s) noexcept { std::swap(d, s.d); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float64; }
        [[nodiscard]] __host__ __device__ double getTrivial() const noexcept { return d; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept{ return d == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return d > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return d < 0; }
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        [[nodiscard]] bool isInteger() const;
        /* Static Members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Scalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
        [[nodiscard]] inline static Scalar random_uniform(RandomPool<RandomGenerator, FixedSeed>& pool);
        template<class RandomGenerator>
        [[nodiscard]] inline static Scalar random_normal(RandomGenerator& gen);
        template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
        [[nodiscard]] inline static Scalar random_normal(RandomPool<RandomGenerator, FixedSeed>& pool);
        template<class RandomPoolType>
        [[nodiscard]] static Scalar random_normal(GaussRandomPool<This, RandomPoolType>& pool) { return pool(); }
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Scalar random_any(Distribution& dist, RandomGenerator& gen);
    #ifdef PHYSICA_HDF5
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_DOUBLE; }
    #endif
    };

    template<class OtherScalar>
    __host__ __device__ inline Scalar<Float64>::Scalar(const ScalarBase<OtherScalar>& s) : d(double(s.getDerived())) {}

    __host__ __device__ inline bool Scalar<Float64>::isFinite() const noexcept {
    #ifdef __CUDA_ARCH__
        return isfinite(d);
    #else
        return std::isfinite(d);
    #endif
    }

    template<class RandomGenerator>
    inline Scalar<Float64> Scalar<Float64>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<double> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline Scalar<Float64> Scalar<Float64>::random_uniform(RandomPool<RandomGenerator, FixedSeed>& pool) {
        return random_uniform<RandomGenerator>(pool.getGen());
    }

    template<class RandomGenerator>
    inline Scalar<Float64> Scalar<Float64>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<double> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline Scalar<Float64> Scalar<Float64>::random_normal(RandomPool<RandomGenerator, FixedSeed>& pool) {
        return random_normal<RandomGenerator>(pool.getGen());
    }

    template<class Distribution, class RandomGenerator>
    inline Scalar<Float64> Scalar<Float64>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Scalar(dist(gen));
    }

    inline std::ostream& operator<<(std::ostream& os, const Scalar<Float64>& s) {
        const auto lastPrec = os.precision();
        return os << std::setprecision(16) << float(s) << std::setprecision(lastPrec);
    }

    inline Scalar<Float64>& operator++(Scalar<Float64>& s) {
        s += float64(1.0);
        return s;
    }

    inline Scalar<Float64>& operator--(Scalar<Float64>& s) {
        s -= float64(1.0);
        return s;
    }

    inline Scalar<Float64> operator++(Scalar<Float64>& s, int) {
        float64 temp(s);
        s += float64(1.0);
        return temp;
    }

    inline Scalar<Float64> operator--(Scalar<Float64>& s, int) {
        float64 temp(s);
        s -= float64(1.0);
        return temp;
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::Float64>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<double>
    #else
        public numeric_limits<double>
    #endif
    {};
}
