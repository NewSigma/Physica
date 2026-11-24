/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
