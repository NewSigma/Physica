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
#pragma once

namespace Physica::Core {
    enum UnitSystem {
        SI,
        AU
    };

    template<UnitSystem unitSystem = SI> class PhyConst;
    /**
     * Physical constants are from [1]
     * Relative atom mass are from internet
     * 
     * Reference:
     * [1] Physical constants, SI (NIST 2018) http://physics.nist.gov/constants
     */
    template<>
    class PhyConst<SI> {
    public:
        constexpr static double planck = 6.62607015E-34;
        constexpr static double reducedPlanck = planck / (2 * M_PI);
        constexpr static double electroMass = 9.1093837015E-31;
        constexpr static double unitCharge = 1.602176634E-19;
        constexpr static double bohrRadius = 5.29177210903E-11;
        constexpr static double protonMass = 1.67262192369E-27;
        constexpr static double neutronMass = 1.67492749804E-27;
        constexpr static double vacuumDielectric = 8.8541878128E-12;
        constexpr static double boltzmannK = 1.380649E-23;
        /**
         * The first element is a space holder
         */
        constexpr static double relativeAtomMass[10]{0, 1.00794, 4.002602, 6.941, 9.012182, 10.806, 12.0096, 14.00643, 15.99903, 18.9984032};
    };

    template<>
    class PhyConst<AU> {
    private:
        constexpr static double hartreeInEv = 27.211386245988;
        constexpr static double rydbergInEv = hartreeInEv * 0.5;
        constexpr static double bohrInAngstorm = PhyConst<SI>::bohrRadius * 1E10;
    public:
        constexpr static double planck = M_PI * 2;
        constexpr static double reducedPlanck = 1;
        constexpr static double electronMass = 1;
        constexpr static double unitCharge = 1;
        constexpr static double bohrRadius = 1;
        constexpr static double protonMass = PhyConst<SI>::protonMass / PhyConst<SI>::electroMass;
        constexpr static double neutronMass = PhyConst<SI>::neutronMass / PhyConst<SI>::electroMass;

        [[nodiscard]] constexpr static double hartreeToEv(double hartree) { return hartree * hartreeInEv; }
        [[nodiscard]] constexpr static double eVToHartree(double ev) { return ev * (1.0 / hartreeInEv); }
        [[nodiscard]] constexpr static double bohrToAngstorm(double bohr) { return bohr * bohrInAngstorm; }
        [[nodiscard]] constexpr static double angstormToBohr(double angstorm) { return angstorm * (1.0 / bohrInAngstorm); }
    };
}