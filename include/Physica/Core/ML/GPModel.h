/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"

namespace Physica {
    /**
     * GP - Gauss Process
     */
    template<Scalar T>
    class GPModel {
        using This = GPModel;
        using Tv = T::ValueType;
    public:
        class GaussKernel;
    private:
        GaussKernel kernel;

        DenseLU<Tv, true> covarLU;
        VectorND<Tv> coeffs;
    public:
        GPModel(size_t numFeature);
        GPModel(const This&) = default;
        GPModel(This&&) noexcept = default;
        ~GPModel() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        CoDiff<T> regression(const MatrixND<Tv>& sampleX, const VectorND<Tv>& sampleY, const auto& uncertainty);
        void step(auto& opt) noexcept { kernel.step(opt); }
        void zero_grad() noexcept { kernel.zero_grad(); }

        [[nodiscard]] Vector2D<Tv> predict(const VectorND<Tv>& x) const;
        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getKernel() const noexcept { return kernel; }
        [[nodiscard]] const auto& getSampleX() const noexcept { return kernel.getSampleX(); }
        [[nodiscard]] size_t getNumFeature() const noexcept { return kernel.getNumFeature(); }
        [[nodiscard]] size_t getNumSamples() const noexcept { return kernel.getOrder(); }
        [[nodiscard]] const auto& getCovarLU() const noexcept { return covarLU; }
    };

    template<Scalar T>
    GPModel<T>::GPModel(size_t numFeature) : kernel(numFeature) {
        assert(numFeature > 0);
    }
    /**
     * \returns marginal likelihood
     */
    template<Scalar T>
    CoDiff<T> GPModel<T>::regression(const MatrixND<Tv>& sampleX, const VectorND<Tv>& sampleY, const auto& uncertainty) {
        assert(sampleX.getRow() == getNumFeature());
        assert(sampleX.getCol() == sampleY.getLength());
        kernel.setSampleX(sampleX);
        covarLU.resize(getNumSamples());

        auto& covars = covarLU.getMatrixLU();
        const Tv var = kernel.getVar();
        assert(!var.isNegative());
        for (size_t major = 0; major < covars.getMaxMajor(); ++major) {
            for (size_t minor = 0; minor < covars.getMaxMinor(); ++minor) {
                bool diag = major == minor;
                if constexpr (Scalar<decltype(uncertainty)>)
                    covars.refFromMajorMinor(major, minor) = diag ? (var + uncertainty) : kernel.calc_value(major, minor);
                else {
                    static_assert(Vector<decltype(uncertainty)>, "[Error]: Unexpected uncertainty type");
                    covars.refFromMajorMinor(major, minor) = diag ? (var + uncertainty.calc(major)) : kernel.calc_value(major, minor);
                }
            }
        }
        covarLU.compute();
        coeffs = covarLU.inv() * sampleY;
        const Tv likelihood = Tv(-0.5) * (sampleY * coeffs + covarLU.lnAbsDet() + ln(Tv(2) * MathConst<Tv>::pi));
        if constexpr (T::isDiffable()) {
            const auto& l = co_yield likelihood;
            kernel.reverse(Tv(2) * l.grad() * reciprocal_elem(coeffs * coeffs.transpose() - MatrixND<Tv>(covarLU.inv())));
        }
        else
            co_return likelihood;
    }

    template<Scalar T>
    auto GPModel<T>::predict(const VectorND<Tv>& x) const -> Vector2D<Tv> {
        assert(x.getLength() == getNumFeature());
        if (getNumSamples() == 0) [[unlikely]]
            return {0, 1};
        auto buffer = VectorND<Tv>::generate([&](size_t i) { return kernel.dot(x, i); }, getNumSamples());
        auto sol = VectorND<Tv>(covarLU.inv() * buffer);
        Tv mean = buffer * coeffs;
        Tv devia = sqrt(std::max(kernel.getVar() - buffer * sol, Tv(0)));
        return {mean, devia};
    }

    template<Scalar T>
    void GPModel<T>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        kernel.swap(obj.kernel);

        covarLU.swap(obj.covarLU);
        coeffs.swap(obj.coeffs);
    }

    template<Scalar T>
    class GPModel<T>::GaussKernel : public RValueMatrix<GaussKernel, T> {
        using This = GaussKernel;
        using Base = RValueMatrix<This, T>;

        MatrixND<Tv> sampleX;
        VectorND<T> alpha;
        T svar;
    public:
        explicit GaussKernel(size_t numFeature);
        GaussKernel(VectorND<Tv> alpha, Tv svar);
        GaussKernel(const This&) = default;
        GaussKernel(This&&) noexcept = default;
        ~GaussKernel() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t r, size_t c) const;
        [[nodiscard]] Tv calc_value(size_t r, size_t c) const;
        [[nodiscard]] Tv dot(const VectorND<Tv>& x, size_t i) const;

        void step(auto& opt) noexcept;
        void zero_grad() noexcept;

        [[nodiscard]] auto&& transpose(this auto&& self) noexcept { return std::forward<decltype(self)>(self); }
        /* Getters */
        [[nodiscard]] const auto& getSampleX() const noexcept { return sampleX; }
        [[nodiscard]] size_t getOrder() const noexcept { return sampleX.getCol(); }
        [[nodiscard]] size_t getNumFeature() const noexcept { return alpha.getLength(); }
        [[nodiscard]] const auto& getAlpha() const noexcept { return alpha; }
        [[nodiscard]] Tv getVar() const noexcept { return abs(svar.value()); }
        /* Setters */
        void setSampleX(const MatrixND<Tv>& sample) { sampleX = sample; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static auto getMajor() noexcept { return MatrixMajor::BothMajor; }
    };

    template<Scalar T>
    GPModel<T>::GaussKernel::GaussKernel(size_t numFeature) : GaussKernel(VectorND<Tv>(numFeature, 1), 1) {}

    template<Scalar T>
    GPModel<T>::GaussKernel::GaussKernel(VectorND<Tv> alpha, Tv svar) : alpha(std::move(alpha)), svar(svar) {
        static_assert(Base::isStaticSquare());
    }

    template<Scalar T>
    CoDiff<T> GPModel<T>::GaussKernel::calc(size_t r, size_t c) const {
        return exp(-square(alpha * (sampleX.col(r) - sampleX.col(c)))) * abs(svar);
    }

    template<Scalar T>
    auto GPModel<T>::GaussKernel::calc_value(size_t r, size_t c) const -> Tv {
        return exp(-square(alpha.values() * (sampleX.col(r) - sampleX.col(c)))) * abs(svar.value());
    }

    template<Scalar T>
    auto GPModel<T>::GaussKernel::dot(const VectorND<Tv>& x, size_t i) const -> Tv {
        return exp(-square(alpha.values() * (x - sampleX.col(i)))) * abs(svar.value());
    }

    template<Scalar T>
    void GPModel<T>::GaussKernel::step(auto& opt) noexcept {
        opt.step(alpha, svar);
    }

    template<Scalar T>
    void GPModel<T>::GaussKernel::zero_grad() noexcept {
        alpha.zero_grad();
        svar.zero_grad();
    }
}
