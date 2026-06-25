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
    public:
        class GaussKernel;
    private:
        GaussKernel kernel;

        DenseLU<T, true> lu;
        MatrixND<T> sampleX;
        VectorND<T> coeffs;
    public:
        GPModel(size_t numFeature);
        GPModel(const This&) = default;
        GPModel(This&&) noexcept = default;
        ~GPModel() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        void regression(const MatrixND<T>& sampleX_, const VectorND<T>& sampleY, const auto& uncertainty);
        [[nodiscard]] Vector2D<T> predict(const VectorND<T>& x) const;

        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumFeature() const noexcept { return kernel.getNumFeature(); }
        [[nodiscard]] size_t getNumSamples() const noexcept { return coeffs.getLength(); }
        [[nodiscard]] const auto& getLU() const noexcept { return lu; }
    };

    template<Scalar T>
    GPModel<T>::GPModel(size_t numFeature) : kernel(numFeature) {
        assert(numFeature > 0);
    }

    template<Scalar T>
    void GPModel<T>::regression(const MatrixND<T>& sampleX_, const VectorND<T>& sampleY, const auto& uncertainty) {
        assert(sampleX_.getRow() == getNumFeature());
        assert(sampleX_.getCol() == sampleY.getLength());
        sampleX = sampleX_;
        lu.resize(sampleX.getCol());

        auto& covars = lu.getMatrixLU();
        for (size_t major = 0; major < covars.getMaxMajor(); ++major) {
            for (size_t minor = 0; minor < covars.getMaxMinor(); ++minor) {
                bool diag = major == minor;
                if constexpr (Scalar<decltype(uncertainty)>)
                    covars.refFromMajorMinor(major, minor) = diag ? (kernel.getVar() + uncertainty) : kernel(sampleX.col(major), sampleX.col(minor));
                else {
                    static_assert(Vector<decltype(uncertainty)>, "[Error]: Unexpected uncertainty type");
                    covars.refFromMajorMinor(major, minor) = diag ? (kernel.getVar() + uncertainty[major]) : kernel(sampleX.col(major), sampleX.col(minor));
                }
            }
        }
        lu.compute();
        coeffs = lu.inv() * sampleY;
    }

    template<Scalar T>
    Vector2D<T> GPModel<T>::predict(const VectorND<T>& x) const {
        assert(x.getLength() == getNumFeature());
        if (getNumSamples() == 0) [[unlikely]]
            return {0, 1};
        auto buffer = VectorND<T>::generate(getNumSamples(), [&](size_t i) { return kernel(x, sampleX.col(i)); });
        auto sol = VectorND<T>(lu.inv() * buffer);
        T mean = buffer * coeffs;
        T devia = sqrt(std::max(kernel.getVar() - buffer * sol, T(0)));
        return {mean, devia};
    }

    template<Scalar T>
    void GPModel<T>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        kernel.swap(obj.kernel);

        lu.swap(obj.lu);
        coeffs.swap(obj.coeffs);
    }

    template<Scalar T>
    class GPModel<T>::GaussKernel {
        using This = GaussKernel;

        VectorND<T> alpha;
        T var;
    public:
        explicit GaussKernel(size_t numFeature);
        GaussKernel(VectorND<T> alpha, T var);
        GaussKernel(const This&) = default;
        GaussKernel(This&&) noexcept = default;
        ~GaussKernel() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] T operator()(const Vector auto& x1, const Vector auto& x2) const;
        /* Getters */
        [[nodiscard]] size_t getNumFeature() const noexcept { return alpha.getLength(); }
        [[nodiscard]] T getVar() const noexcept { return var; }
    };

    template<Scalar T>
    GPModel<T>::GaussKernel::GaussKernel(size_t numFeature) : GaussKernel(VectorND<T>(numFeature, 1), 1) {}

    template<Scalar T>
    GPModel<T>::GaussKernel::GaussKernel(VectorND<T> alpha, T var) : alpha(std::move(alpha)), var(var) {}

    template<Scalar T>
    T GPModel<T>::GaussKernel::operator()(const Vector auto& x1, const Vector auto& x2) const {
        return var * exp(-alpha * square(x1 - x2));
    }
}
