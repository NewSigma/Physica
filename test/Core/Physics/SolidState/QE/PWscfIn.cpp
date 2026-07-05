/*
 * Copyright 2023-2025 Weibo He.
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
#include <fstream>
#include "Physica/Core/Physics/SolidState/QE/PWscfIn.h"
#include "Physica/Core/Utils/Unix/TempFile.h"

using namespace Physica;

namespace {
    const char* data = "&CONTROL\n"
                    "  pseudo_dir = '/home/sigma/Program/pot-QE/pslibrary-master/pbe/PSEUDOPOTENTIALS'\n"
                    "  calculation = 'scf'\n"
                    "  prefix = 'fBH'\n"
                    "  outdir = './output'\n"
                    "  tprnfor = .true.\n"
                    "  tstress = .true.\n"
                    "/\n"
                    "\n"
                    "&SYSTEM\n"
                    "  nat = 3\n"
                    "  ntyp = 2\n"
                    "  ibrav = 0\n"
                    "  ecutwfc = 60\n"
                    "  ecutrho = 480\n"
                    "  occupations = 'fixed'\n"
                    "/\n"
                    "\n"
                    "&ELECTRONS\n"
                    "  mixing_beta = 0.7\n"
                    "  conv_thr =  1.0d-8\n"
                    "/\n"
                    "\n"
                    "ATOMIC_SPECIES\n"
                    "H   1.00798  H.pbe-kjpaw_psl.1.0.0.UPF\n"
                    "O   15.9994  O.pbe-nl-kjpaw_psl.1.0.0.UPF\n"
                    "\n"
                    "CELL_PARAMETERS (angstrom)\n"
                    "   4.648733580  -0.001573654   0.000000000\n"
                    "   2.326056755   4.126723845   0.000000000\n"
                    "   0.000000000   0.000000000  18.000000000\n"
                    "\n"
                    "ATOMIC_POSITIONS (angstrom)\n"
                    "H             3.8682839689        1.6444713403        3.1314180072\n"
                    "H             3.8678304286        0.6390600192        1.9183297424\n"
                    "O             3.8678914803        0.6699600914        2.9195167806\n"
                    "\n"
                    "K_POINTS automatic\n"
                    "4 4 1 0 0 0\n";
}


int main() {
    try {
        auto tmp = TempFile("/tmp/tmpXXXXXX");
        {
            std::ofstream fout(tmp.getName());
            fout << data;
        }
        std::ifstream fin(tmp.getName());
        PWscfIn input{};
        fin >> input;
    }
    catch (std::exception& e) {
        std::cout << e.what() << '\n';
        return 1;
    }
    return 0;
}
