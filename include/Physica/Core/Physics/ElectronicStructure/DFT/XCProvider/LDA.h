/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Physics/ElectronicStructure/DFT/Grid3D.h"

namespace Physica::Core {
    enum class LDAType {
        HL
    };

    template<class ScalarType, LDAType type, bool isSpinPolarized> class LDA;

    namespace Internal {
        template<class T> class Traits;

        template<class ScalarType, LDAType type, bool polarized>
        class Traits<LDA<ScalarType, type, polarized>> {
            using Grid = Grid3D<ScalarType, false>;
        public:
            constexpr static bool isSpinPolarized = polarized;
            using DensityType = typename std::conditional<isSpinPolarized, std::pair<Grid, Grid>, Grid>::type;
            using PotType = typename std::conditional<isSpinPolarized, std::pair<Grid, Grid>, Grid>::type;
        };
    }

    template<class ScalarType, LDAType type>
    class LDA<ScalarType, type, true> {
        using This = LDA<ScalarType, type, true>;
    public:
        constexpr static bool isSpinPolarized = Internal::Traits<This>::isSpinPolarized;
        using DensityType = typename Internal::Traits<This>::DensityType;
        using PotType = typename Internal::Traits<This>::PotType;
    private:
        Vector<ScalarType> buffer;
        Vector<ScalarType> buffer1;
        Vector<ScalarType> buffer2;
        Vector<ScalarType> buffer3;
    public:
        LDA(size_t bufferSize);
        LDA(const LDA&) = delete;
        LDA(LDA&&) = delete;
        ~LDA() = default;
        /* Operators */
        LDA& operator=(const LDA&) = delete;
        LDA& operator=(LDA&&) = delete;
        /* Operations */
        void fill(const DensityType& density, PotType& xc);
        /* Getters */
        [[nodiscard]] size_t getBufferSize() const noexcept { return buffer.getLength(); }
    private:
        void fillExchange(const DensityType& density, PotType& xc);
        void addCorreclation(const DensityType& density, PotType& xc);
    };

    template<class ScalarType, LDAType type>
    LDA<ScalarType, type, true>::LDA(size_t bufferSize) : buffer(bufferSize)
                                                        , buffer1(bufferSize)
                                                        , buffer2(bufferSize)
                                                        , buffer3(bufferSize) {}

    template<class ScalarType, LDAType type>
    void LDA<ScalarType, type, true>::fill(const DensityType& density, PotType& xc) {
        [[maybe_unused]] const auto& rho = density.first.asVector();
        [[maybe_unused]] const auto& zeta = density.second.asVector();
        assert(rho.getLength() == getBufferSize());
        assert(zeta.getLength() == getBufferSize());
        fillExchange(density, xc);
        addCorreclation(density, xc);
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:106
     */
    template<class ScalarType, LDAType type>
    void LDA<ScalarType, type, true>::fillExchange(const DensityType& density, PotType& xc) {
        constexpr double factor0 = -0.73855876638202240588;
        constexpr double factor1 = -0.93052573634910002500;
        constexpr double factor_f = 1.9236610509315363198;
        constexpr double factor_df = -2.5648814012420484263;
        const auto& rho = density.first.asVector();
        const auto& zeta = density.second.asVector();
        auto& xc_up = xc.first.asVector();
        auto& xc_down = xc.second.asVector();
        auto& buffer5 = xc_up;
        auto& buffer6 = xc_down;

        buffer = pow(rho, ScalarType(1.0 / 3));
        buffer1 = -zeta + ScalarType::One();
        buffer2 = zeta + ScalarType::One();
        buffer3 = pow(buffer1, ScalarType(1.0 / 3));
        buffer5 = hadamard(buffer3, buffer3);
        buffer5 = hadamard(buffer5, buffer5);

        buffer6 = pow(buffer2, ScalarType(1.0 / 3));
        auto& df = buffer3;
        df = (buffer3 - buffer6) * ScalarType(factor_df);
        buffer6 = hadamard(buffer6, buffer6);
        buffer6 = hadamard(buffer6, buffer6);

        auto& unnormalized_f = buffer5;
        unnormalized_f = buffer5 + buffer6 - ScalarType::Two();
        auto& epsilon = buffer5;
        epsilon = ScalarType(4.0 / 3) * hadamard(ScalarType(factor0) + unnormalized_f * ScalarType(factor_f * (factor1 - factor0)), buffer);
        buffer6 = epsilon;

        xc_up += hadamard(hadamard(buffer1, df), buffer) * ScalarType(factor1 - factor0);
        xc_down += hadamard(hadamard(buffer2, df), buffer) * ScalarType(factor0 - factor1);
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:479-480
     */
    template<class ScalarType, LDAType type>
    void LDA<ScalarType, type, true>::addCorreclation([[maybe_unused]] const DensityType& density, PotType& xc) {
        auto& xc_up = xc.first.asVector();
        auto& xc_down = xc.second.asVector();
        switch(type) {
        case LDAType::HL:
            constexpr double factor1 = -0.045 / 2;
            constexpr double factor2 = 33.851831034345862;
            buffer = ScalarType(factor1) * ln(ScalarType::One() + ScalarType(factor2) * buffer);
            xc_up += buffer;
            xc_down += buffer;
        }
    }

    template<class ScalarType, LDAType type>
    class LDA<ScalarType, type, false> {
        using This = LDA<ScalarType, type, false>;
    public:
        constexpr static bool isSpinPolarized = Internal::Traits<This>::isSpinPolarized;
        using DensityType = typename Internal::Traits<This>::DensityType;
        using PotType = typename Internal::Traits<This>::PotType;
    private:
        Vector<ScalarType> buffer;
    public:
        LDA(size_t bufferSize);
        LDA(const LDA&) = delete;
        LDA(LDA&&) = delete;
        ~LDA() = default;
        /* Operators */
        LDA& operator=(const LDA&) = delete;
        LDA& operator=(LDA&&) = delete;
        /* Operations */
        void fill(const DensityType& density, PotType& xc);
        /* Getters */
        [[nodiscard]] size_t getBufferSize() const noexcept { return buffer.getLength(); }
    private:
        void fillExchange(const DensityType& density, PotType& xc);
        void addCorreclation(const DensityType& density, PotType& xc);
    };

    template<class ScalarType, LDAType type>
    LDA<ScalarType, type, false>::LDA(size_t bufferSize) : buffer(bufferSize) {}

    template<class ScalarType, LDAType type>
    void LDA<ScalarType, type, false>::fill(const DensityType& density, PotType& xc) {
        [[maybe_unused]] const auto& rho = density.asVector();
        assert(rho.getLength() == getBufferSize());
        fillExchange(density, xc);
        addCorreclation(density, xc);
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:106
     */
    template<class ScalarType, LDAType type>
    void LDA<ScalarType, type, false>::fillExchange(const DensityType& density, PotType& xc) {
        constexpr double factor0 = -0.73855876638202240588;
        const auto& rho = density.asVector();
        auto& V = xc.asVector();

        buffer = pow(rho, ScalarType(1.0 / 3));
        V = ScalarType(4.0 / 3 * factor0) * buffer;
    }
    /**
     * Reference:
     * [1] Martin,Richard M. Electronic structure : basic theory and practical methods[M].Beijing: World publishing corporation; Cambridge: Cambridge University Press, 2017:479-480
     */
    template<class ScalarType, LDAType type>
    void LDA<ScalarType, type, false>::addCorreclation([[maybe_unused]] const DensityType& density, PotType& xc) {
        auto& V = xc.asVector();
        switch(type) {
        case LDAType::HL:
            constexpr double factor1 = -0.045 / 2;
            constexpr double factor2 = 33.851831034345862;
            buffer = ScalarType(factor1) * ln(ScalarType::One() + ScalarType(factor2) * buffer);
            V += buffer;
        }
    }
}
