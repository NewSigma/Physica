/*
 * Copyright 2023 WeiBo He.
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
#include <iostream>
#include "Physica/Core/IO/QE/PWscfIn.h"
#include "Physica/Core/Exception/BadFileFormatException.h"

namespace Physica::Core {
    PWscfIn& PWscfIn::operator=(PWscfIn obj) noexcept {
        swap(obj);
        return *this;
    }

    std::ostream& operator<<(std::ostream& os, const PWscfIn& input) {
        os << "&CONTROL\n";
        {
            os << "  calculation = '" << PWscfIn::calculationToStr(input.calculation) << "'\n";
            os << "  tstress = ." << (input.tstress ? "true" : "false") << ".\n";
            os << "  tprnfor = ." << (input.tprnfor ? "true" : "false") << ".\n";
            os << "  outdir = '" << input.outdir.c_str() << "'\n";
            os << "  prefix = '" << input.prefix.c_str() << "'\n";
            os << "  pseudo_dir = '" << input.pseudo_dir.c_str() << "'\n";
        }
        os << "/\n";
        return os;
    }

    std::istream& operator>>(std::istream& is, PWscfIn& input) {
        assert(is.good());
        Utils::Array<char> buffer(PWscfIn::BufferSize);
        do {
            is.getline(buffer.data(), buffer.size());
            const std::string card = buffer.data();
            const bool isControl = card == "&CONTROL";
            if (isControl)
                input.readControl(is, buffer);
        } while(is.peek() != EOF);

        if (!is)
            throw BadFileFormatException("[Error]: Bad PWscf input file");
        return is;
    }

    void PWscfIn::swap(PWscfIn& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(calculation, obj.calculation);
        std::swap(tstress, obj.tstress);
        std::swap(tprnfor, obj.tprnfor);
        outdir.swap(obj.outdir);
        prefix.swap(obj.prefix);
        pseudo_dir.swap(obj.pseudo_dir);
    }

    const char* PWscfIn::calculationToStr(CalculationType calculation) {
        switch (calculation) {
        case SCF:
            return "scf";
        case NSCF:
            return "nscf";
        case Bands:
            return "bands";
        case Relax:
            return "relax";
        case MD:
            return "md";
        case VC_Relax:
            return "vc-relax";
        case VC_MD:
            return "vc_md";
        default:
            assert(calculation == SCF && "[Error]: Invalid calculation type");
            return "[Error]: Invalid calculation type";
        };
    }

    void PWscfIn::readControl(std::istream& is, Utils::Array<char>& buffer) {
        do {
            is >> std::ws;
            is.getline(buffer.data(), buffer.getLength(), ' ');
            const std::string tag = buffer.data();
            if (tag == "calculation") {
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\'');
                is.getline(buffer.data(), buffer.size(), '\'');
                setCalculation(buffer.data());
            }
            else if (tag == "tstress")
                tstress = readBool(is, buffer);
            else if (tag == "tprnfor")
                tprnfor = readBool(is, buffer);
            else if (tag == "outdir")
                readStr(is, buffer, outdir);
            else if (tag == "prefix")
                readStr(is, buffer, prefix);
            else if (tag == "pseudo_dir")
                readStr(is, buffer, pseudo_dir);
            else
                throw BadFileFormatException("[Error]: Bad tag in CONTROL");
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } while(is.peek() != '/');
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    void PWscfIn::setCalculation(const std::string& str) {
        if (str == "scf")
            calculation = SCF;
        else if (str == "nscf")
            calculation = NSCF;
        else if (str == "bands")
            calculation = Bands;
        else if (str == "relax")
            calculation = Relax;
        else if (str == "md")
            calculation = MD;
        else if (str == "vc-relax")
            calculation = VC_Relax;
        else if (str == "vc-md")
            calculation = VC_MD;
    }

    void PWscfIn::readStr(std::istream& is, Utils::Array<char>& buffer, std::string& saveTo) {
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\'');
        is.getline(buffer.data(), buffer.getLength(), '\'');
        saveTo = buffer.data();
    }

    bool PWscfIn::readBool(std::istream& is, Utils::Array<char>& buffer) {
        is.ignore(std::numeric_limits<std::streamsize>::max(), '.');
        is.getline(buffer.data(), buffer.getLength(), '.');
        return strcmp(buffer.data(), "true") == 0;
    }
}
