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

namespace Physica {
    template<>
    class Traits<Core::float32> {
    public:
        using ScalarType = Core::float32;
        using RealType = ScalarType;
        using ComplexType = Core::ComplexScalar<ScalarType>;
        using TrivialType = float;
        using PlainScalar = ScalarType;
        constexpr static Core::ScalarOption Option = Core::Float32;
        constexpr static bool isComplex = false;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;
        constexpr static unsigned int Order = 0;
    };
}

namespace Physica::Core {
    template<>
    class Scalar<Float> : public ScalarBase<Scalar<Float>> {
        using This = Scalar<Float>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        float f;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(float f_) : f(f_) {}
        Scalar(const Integer& i) : Scalar(float(double(i))) {}
        Scalar(const Rational& r) : Scalar(float(double(r))) {}
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
        __host__ __device__ explicit operator float() const { return f; }
        __host__ __device__ explicit operator double() const { return f; }
        __host__ __device__ Scalar operator+(const Scalar& s) const { return Scalar(f + s.f); }
        __host__ __device__ Scalar operator-(const Scalar& s) const { return Scalar(f - s.f); }
        __host__ __device__ Scalar operator*(const Scalar& s) const { return Scalar(f * s.f); }
        __host__ __device__ Scalar operator/(const Scalar& s) const { return Scalar(f / s.f); }
        Scalar operator<<(int i) const { return Scalar(f * std::pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(f / std::pow(2, i)); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-f); }
        __host__ __device__ bool operator>(const Scalar& s) const { return f > s.f; }
        __host__ __device__ bool operator<(const Scalar& s) const { return f < s.f; }
        __host__ __device__ bool operator==(const Scalar& s) const { return f == s.f; }
        PHYSICA_API friend std::istream& operator>>(std::istream& is, Scalar& scalar);
        /* Operations */
        void swap(Scalar& __restrict s) noexcept { std::swap(f, s.f); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float; }
        [[nodiscard]] __host__ __device__ float getTrivial() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return f == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return f > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return f < 0; }
        [[nodiscard]] bool isFinite() const noexcept { return std::isfinite(f); }
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
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_FLOAT; }
    #endif
    };

    template<class OtherScalar>
    __host__ __device__ inline Scalar<Float>::Scalar(const ScalarBase<OtherScalar>& s) : f(float(s.getDerived())) {}

    template<class RandomGenerator>
    inline Scalar<Float> Scalar<Float>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<float> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline Scalar<Float> Scalar<Float>::random_uniform(RandomPool<RandomGenerator, FixedSeed>& pool) {
        return random_uniform<RandomGenerator>(pool.getGen());
    }

    template<class RandomGenerator>
    inline Scalar<Float> Scalar<Float>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<float> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline Scalar<Float> Scalar<Float>::random_normal(RandomPool<RandomGenerator, FixedSeed>& pool) {
        return random_normal<RandomGenerator>(pool.getGen());
    }

    template<class Distribution, class RandomGenerator>
    inline Scalar<Float> Scalar<Float>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Scalar(dist(gen));
    }

    inline std::ostream& operator<<(std::ostream& os, const Scalar<Float32>& s) {
        const auto lastPrec = os.precision();
        return os << std::setprecision(7) << float(s) << std::setprecision(lastPrec);
    }

    inline Scalar<Float>& operator++(Scalar<Float>& s) {
        s += Scalar<Float>(1.0F);
        return s;
    }

    inline Scalar<Float>& operator--(Scalar<Float>& s) {
        s -= Scalar<Float>(1.0F);
        return s;
    }

    inline Scalar<Float> operator++(Scalar<Float>& s, int) {
        Scalar<Float> temp(s);
        s += Scalar<Float>(1.0F);
        return temp;
    }

    inline Scalar<Float> operator--(Scalar<Float>& s, int) {
        Scalar<Float> temp(s);
        s -= Scalar<Float>(1.0F);
        return temp;
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::Float>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<float>
    #else
        public numeric_limits<float>
    #endif
    {};
}
