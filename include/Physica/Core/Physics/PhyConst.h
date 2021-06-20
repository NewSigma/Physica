/*
 * Copyright 2021 WeiBo He.
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
namespace Physica::Core::Physics {
    /**
     * Physical constants in International System of Units
     */
    class PhyConst {
    public:
        const static double planck = 6.6260755E-34;
        const static double reducedPlanck = 1.05457266E-34;
        const static double electroMass = 9.10938215E-31;
        const static double protonMass = 1.6726231E-27;
        const static double neutronMass = 1.6749286E-27;
    };
}