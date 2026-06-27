/*
 * Copyright 2023-2026 Weibo He.
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

#include <iosfwd>
#include <string>
#include "Physica/Core/Scalar/Real.h"

namespace Physica {
    /**
     * Reference:
     * [1] https://www.quantum-espresso.org/Doc/INPUT_PW.html
     */
    class PHYSICA_API PWscfIn {
        using This = PWscfIn;
        using ScalarType = float64;
        constexpr static size_t BufferSize = 1024;
    public:
        enum CalculationType {
            SCF,
            NSCF,
            Bands,
            Relax,
            MD,
            VC_Relax,
            VC_MD
        };
    private:
        CalculationType calculation;
        bool tstress;
        bool tprnfor;
        Array<char> outdir;
        Array<char> prefix;
        Array<char> pseudo_dir;
    public:
        PWscfIn() = default;
        PWscfIn(const This&) = default;
        PWscfIn(This&&) noexcept = default;
        ~PWscfIn() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        friend std::ostream& operator<<(std::ostream& os, const This& input);
        friend std::istream& operator>>(std::istream& is, This& input);
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        static const char* calculationToStr(CalculationType calculation);
    private:
        void readControl(std::istream& is, Array<char>& buffer);
        void setCalculation(const std::string& str);
        void readStr(std::istream& is, Array<char>& buffer, Array<char>& saveTo);
        bool readBool(std::istream& is, Array<char>& buffer);
    };

    PHYSICA_API std::ostream& operator<<(std::ostream& os, const PWscfIn& input);
    PHYSICA_API std::istream& operator>>(std::istream& is, PWscfIn& input);
}
