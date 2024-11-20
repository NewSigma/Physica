/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h"

namespace Physica::Core {
    enum class LDAType {
        HL
    };

    template<Scalar T, LDAType type, bool IsSpinPolarized> class LDA;

    template<Scalar T, LDAType type>
    class LDA<T, type, true> {
        using This = LDA<T, type, true>;
    public:
        constexpr static bool IsSpinPolarized = Traits<This>::IsSpinPolarized;
        using DensityType = typename Traits<This>::DensityType;
        using PotType = typename Traits<This>::PotType;
    private:
        VectorND<T> buffer;
        VectorND<T> buffer1;
        VectorND<T> buffer2;
        VectorND<T> buffer3;
    public:
        LDA() = default;
        LDA(size_t bufferSize);
        LDA(const LDA&) = default;
        LDA(LDA&&) = default;
        ~LDA() = default;
        /* Operators */
        LDA& operator=(LDA lda) noexcept;
        /* Operations */
        void fill(const DensityType& density, PotType& xc);
        void swap(LDA& __restrict lda) noexcept;
        /* Getters */
        [[nodiscard]] size_t getBufferSize() const noexcept { return buffer.getLength(); }
    private:
        void fillExchange(const DensityType& density, PotType& xc);
        void addCorreclation(const DensityType& density, PotType& xc);
    };

    template<Scalar T, LDAType type>
    LDA<T, type, true>::LDA(size_t bufferSize) : buffer(bufferSize)
                                                        , buffer1(bufferSize)
                                                        , buffer2(bufferSize)
                                                        , buffer3(bufferSize) {}

    template<Scalar T, LDAType type>
    LDA<T, type, true>& LDA<T, type, true>::operator=(LDA lda) noexcept {
        swap(lda);
        return *this;
    }

    template<Scalar T, LDAType type>
    void LDA<T, type, true>::fill(const DensityType& density, PotType& xc) {
        [[maybe_unused]] const auto& rho = density[0].flatten();
        [[maybe_unused]] const auto& zeta = density[1].flatten();
        assert(rho.getLength() == getBufferSize());
        assert(zeta.getLength() == getBufferSize());
        fillExchange(density, xc);
        addCorreclation(density, xc);
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:106
     */
    template<Scalar T, LDAType type>
    void LDA<T, type, true>::fillExchange(const DensityType& density, PotType& xc) {
        constexpr double factor0 = -0.73855876638202240588;
        constexpr double factor1 = -0.93052573634910002500;
        constexpr double factor_f = 1.9236610509315363198;
        constexpr double factor_df = -2.5648814012420484263;
        const auto& rho = density[0].flatten();
        const auto& zeta = density[1].flatten();
        auto xc_up = xc[0].flatten();
        auto xc_down = xc[1].flatten();
        auto& buffer5 = xc_up;
        auto& buffer6 = xc_down;

        buffer = pow(rho, T(1.0 / 3));
        buffer1 = -zeta + T(1);
        buffer2 = zeta + T(1);
        buffer3 = pow(buffer1, T(1.0 / 3));
        buffer5 = hadamard(buffer3, buffer3);
        buffer5 = hadamard(buffer5, buffer5);

        buffer6 = pow(buffer2, T(1.0 / 3));
        auto& df = buffer3;
        df = (buffer3 - buffer6) * T(factor_df);
        buffer6 = hadamard(buffer6, buffer6);
        buffer6 = hadamard(buffer6, buffer6);

        auto& unnormalized_f = buffer5;
        unnormalized_f = buffer5 + buffer6 - T(2);
        auto& epsilon = buffer5;
        epsilon = T(4.0 / 3) * hadamard(T(factor0) + unnormalized_f * T(factor_f * (factor1 - factor0)), buffer);
        buffer6 = epsilon;

        xc_up += hadamard(hadamard(buffer1, df), buffer) * T(factor1 - factor0);
        xc_down += hadamard(hadamard(buffer2, df), buffer) * T(factor0 - factor1);
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:479-480
     */
    template<Scalar T, LDAType type>
    void LDA<T, type, true>::addCorreclation([[maybe_unused]] const DensityType& density, PotType& xc) {
        auto xc_up = xc[0].flatten();
        auto xc_down = xc[1].flatten();
        switch(type) {
        case LDAType::HL:
            constexpr double factor1 = -0.045 / 2;
            constexpr double factor2 = 33.851831034345862;
            buffer = T(factor1) * ln(T(1) + T(factor2) * buffer);
            xc_up += buffer;
            xc_down += buffer;
        }
    }

    template<Scalar T, LDAType type>
    void LDA<T, type, true>::swap(LDA& __restrict lda) noexcept {
        assert(this != &lda && "[Error]: Self swap is likely a bug");
        buffer.swap(lda.buffer);
        buffer1.swap(lda.buffer1);
        buffer2.swap(lda.buffer2);
        buffer3.swap(lda.buffer3);
    }

    template<Scalar T, LDAType type>
    class LDA<T, type, false> {
        using This = LDA<T, type, false>;
    public:
        constexpr static bool IsSpinPolarized = Traits<This>::IsSpinPolarized;
        using DensityType = typename Traits<This>::DensityType;
        using PotType = typename Traits<This>::PotType;
    private:
        VectorND<T> buffer;
    public:
        LDA() = default;
        LDA(size_t bufferSize);
        LDA(const LDA&) = default;
        LDA(LDA&&) = default;
        ~LDA() = default;
        /* Operators */
        LDA& operator=(LDA lda) noexcept;
        /* Operations */
        void fill(const DensityType& density, PotType& xc);
        void swap(LDA& __restrict lda) noexcept;
        /* Getters */
        [[nodiscard]] size_t getBufferSize() const noexcept { return buffer.getLength(); }
    private:
        void fillExchange(const DensityType& density, PotType& xc);
        void addCorreclation(const DensityType& density, PotType& xc);
    };

    template<Scalar T, LDAType type>
    LDA<T, type, false>::LDA(size_t bufferSize) : buffer(bufferSize) {}

    template<Scalar T, LDAType type>
    LDA<T, type, false>& LDA<T, type, false>::operator=(LDA lda) noexcept {
        swap(lda);
        return *this;
    }

    template<Scalar T, LDAType type>
    void LDA<T, type, false>::fill(const DensityType& density, PotType& xc) {
        [[maybe_unused]] const auto& rho = density[0].flatten();
        assert(rho.getLength() == getBufferSize());
        fillExchange(density, xc);
        addCorreclation(density, xc);
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:106
     */
    template<Scalar T, LDAType type>
    void LDA<T, type, false>::fillExchange(const DensityType& density, PotType& xc) {
        constexpr double factor0 = -0.73855876638202240588;
        const auto& rho = density[0].flatten();
        auto V = xc[0].flatten();

        buffer = cbrt(rho);
        V = T(4.0 / 3 * factor0) * buffer;
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:479-480
     */
    template<Scalar T, LDAType type>
    void LDA<T, type, false>::addCorreclation([[maybe_unused]] const DensityType& density, PotType& xc) {
        auto V = xc[0].flatten();
        switch(type) {
        case LDAType::HL:
            constexpr double factor1 = -0.045 / 2;
            constexpr double factor2 = 33.851831034345862;
            buffer = T(factor1) * ln(T(1) + T(factor2) * buffer);
            V += buffer;
        }
    }

    template<Scalar T, LDAType type>
    void LDA<T, type, false>::swap(LDA& __restrict lda) noexcept {
        assert(this != &lda && "[Error]: Self swap is likely a bug");
        buffer.swap(lda.buffer);
    }
}

namespace Physica {
    template<Scalar T, LDAType type, bool polarized>
    class Traits<LDA<T, type, polarized>> {
        using GridType = RSpaceGrid<T>;
    public:
        constexpr static bool IsSpinPolarized = polarized;
        constexpr static size_t NumSpin = IsSpinPolarized ? 2 : 1;
        using DensityType = Array<GridType, NumSpin>;
        using PotType = Array<GridType, NumSpin>;
    };
}
