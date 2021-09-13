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
        constexpr static double planck = 6.6260755E-34;
        constexpr static double reducedPlanck = 1.05457266E-34;
        constexpr static double electroMass = 9.10938215E-31;
        constexpr static double unitCharge = 1.602176634E-19; 
        constexpr static double protonMass = 1.6726231E-27;
        constexpr static double neutronMass = 1.6749286E-27;
        /**
         * The first element is a space holder
         */
        constexpr static double relativeAtomMass[10]{0, 1.00794, 4.002602, 6.941, 9.012182, 10.806, 12.0096, 14.00643, 15.99903, 18.9984032};
    };
}