/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/Physics/SolidState/IceGenerator.h"

using namespace Physica;
using namespace Physica;
using ScalarType = float64;

namespace Physica {
    class Test {
    public:
        static auto findBondedH(const IceGenerator<ScalarType>& gen, size_t indexO) { return gen.findBondedH(indexO); }
    };
}

CrystalCell<ScalarType> makeCell() {
    using LatticeMatrix = CrystalCell<ScalarType>::LatticeMatrix;
    using PositionMatrix = CrystalCell<ScalarType>::PositionMatrix;
    const LatticeMatrix lattice{
        -5.1399889103954424,    0.0043684262340578,    0.0000000000000000,
        -0.0032117921883728,   -4.1600557380116499,    0.0000000000000000,
        -3.0735325252647820,   -0.7268177903616805,   18.0000000000000000
    };
    const PositionMatrix pos {
        0.0726081710567296, 0.3720055028581801, 0.1023778332101003,
        0.7261072293654427, 0.7598698062359442, 0.0249099824708654,
        0.0780796208184085, 0.0887316261990637,-0.0010676666526211,
        0.8014375091805590, 0.2588755409266756, 0.0305438697861024,
        0.4656459841168016, 0.2719702241499968, 0.9474437306732401,
        0.4548704645367646, 0.8981120521536182, 0.9530859820000415,
        0.4494699053814653, 0.5786334574671179, 0.0565389960871029,
        0.0619486126786336, 0.7438923669781803, 0.1080560466235230,
        0.1704575124358020, 0.5568081344822839, 0.1166696814548458,
        0.9205753978024702, 0.0770597907999365, 0.0405344961066411,
        0.6069508593787584, 0.5815283677184969, 0.0149318696265062,
        0.3570788249281732, 0.0879730900605094, 0.9388101767036668
    };

    CrystalCell<ScalarType> cell({lattice, pos, CrystalCell<ScalarType>::Type::Direct}, {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8});
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    cell.toCartesian();
    return cell;
}

int main() {
    IceGenerator<ScalarType> gen(makeCell(), 5.555, 2.08);
    const auto arr = Test::findBondedH(gen, 0);
    const bool b1 = std::find(arr.cbegin(), arr.cend(), 0) != arr.cend();
    const bool b2 = std::find(arr.cbegin(), arr.cend(), 6) != arr.cend();
    const bool b3 = std::find(arr.cbegin(), arr.cend(), 7) != arr.cend();
    if (!(b1 && b2 && b3))
        return 1;
    return 0;
}
