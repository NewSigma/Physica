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

#include <cmath>

namespace Physica::Core {
    enum UnitSystem {
        SI,
        AU,
        ESU
    };

    template<UnitSystem unitSystem = SI> class PhyConst;
    /**
     * Physical constants are from [1]
     * Relative atom mass are from [2]
     * 
     * Note:
     * Relative atom mass are taken to be the average value of isotopes that are common in nature as [2] does. If you are interesed in
     * a particular isotope, modifications are necessary.
     * 
     * Reference:
     * [1] Physical constants, SI (NIST 2018) http://physics.nist.gov/constants
     * [2] Standard Atomic Weights 2021 https://www.ciaaw.org/atomic-weights.htm
     */
    template<>
    class PhyConst<SI> {
        constexpr static double calorieInJoule = 4.184; // Reference: http://simple.wikipedia.org/wiki/Calorie
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
        constexpr static double avogadroNa = 6.02214076E23;
        constexpr static double speedOfLight = 299792458;
        /**
         * The first element is a space holder
         */
        constexpr static double relativeMassInKg = 1E-3 / avogadroNa;
        constexpr static double relativeAtomMass[19]{0, 1.00798, 4.002602,
                6.968, 9.0121831, 10.814, 12.0106, 14.00686, 15.99940, 18.998403162, 20.1797,
                22.98976928, 24.306, 26.9815384, 28.085, 30.973761998, 32.068, 35.452, 39.878};
        constexpr static const char* elementSymbol[21]{"-", "H", "He",
                "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                "K", "Ca"};
    
        [[nodiscard]] constexpr static double degreeToRadian(double degree) { return degree / 180 * M_PI; }
        [[nodiscard]] constexpr static double atomMass(size_t atomicNum) { return relativeAtomMass[atomicNum] * relativeMassInKg; }
        [[nodiscard]] constexpr static double calorieToJoule(double calorie) { return calorie * calorieInJoule; }
    };
    /**
     * Hartree atomic units
     */
    template<>
    class PhyConst<AU> {
    private:
        constexpr static double hartreeInEv = 27.211386245988;
        constexpr static double rydbergInEv = hartreeInEv * 0.5;
        constexpr static double bohrInAngstorm = PhyConst<SI>::bohrRadius * 1E10;
        constexpr static double hartreeInJoule = PhyConst<SI>::unitCharge * hartreeInEv;
        constexpr static double jouleInHartree = 1 / hartreeInJoule;
        constexpr static double timeInSecond = PhyConst<SI>::reducedPlanck * jouleInHartree;
        constexpr static double temperatureInK = 1 / (PhyConst<SI>::boltzmannK * jouleInHartree);
        constexpr static double pressInGPa = hartreeInJoule / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius) * 1E-9;
    public:
        constexpr static double planck = M_PI * 2;
        constexpr static double reducedPlanck = 1;
        constexpr static double electronMass = 1;
        constexpr static double unitCharge = 1;
        constexpr static double bohrRadius = 1;
        constexpr static double protonMass = PhyConst<SI>::protonMass / PhyConst<SI>::electroMass;
        constexpr static double neutronMass = PhyConst<SI>::neutronMass / PhyConst<SI>::electroMass;
        constexpr static double vacuumDielectric = (PhyConst<SI>::vacuumDielectric * PhyConst<SI>::bohrRadius) / (jouleInHartree * PhyConst<SI>::unitCharge * PhyConst<SI>::unitCharge);
        constexpr static double boltzmannK = 1;
        constexpr static double speedOfLight = PhyConst<SI>::speedOfLight / PhyConst<SI>::bohrRadius * timeInSecond;

        [[nodiscard]] constexpr static double hartreeToEv(double hartree) { return hartree * hartreeInEv; }
        [[nodiscard]] constexpr static double eVToHartree(double ev) { return ev * (1.0 / hartreeInEv); }
        [[nodiscard]] constexpr static double bohrToAngstorm(double bohr) { return bohr * bohrInAngstorm; }
        [[nodiscard]] constexpr static double angstormToBohr(double angstorm) { return angstorm * (1.0 / bohrInAngstorm); }
        [[nodiscard]] constexpr static double timeToSecond(double atomic_time) { return atomic_time * timeInSecond; }
        [[nodiscard]] constexpr static double secondToTime(double second) { return second / timeInSecond; }
        [[nodiscard]] constexpr static double temperatureToK(double atomic_tem) { return atomic_tem * temperatureInK; }
        [[nodiscard]] constexpr static double kToTemperature(double kelvin) { return kelvin / temperatureInK; }
        [[nodiscard]] constexpr static double pressToGPa(double atomic_press) { return atomic_press * pressInGPa; }
        [[nodiscard]] constexpr static double atomMass(size_t atomicNum) { return PhyConst<SI>::atomMass(atomicNum) / PhyConst<SI>::electroMass; }
    };

    template<>
    class PhyConst<ESU> {
        constexpr static double esuToCoulombFactor = 3.335640951084828E-10; // std::sqrt(PhyConst<SI>::vacuumDielectric * 4 * M_PI * 10E-9);
    public:
        constexpr static double unitCharge = PhyConst<SI>::unitCharge / esuToCoulombFactor;
    private:
        constexpr static double debyeToDipoleAUFactor = (10E-10 / unitCharge) * PhyConst<AU>::angstormToBohr(1);
    public:
        [[nodiscard]] constexpr static double esuToCoulomb(double esu) { return esu * esuToCoulombFactor; }
        [[nodiscard]] constexpr static double debyeToDipoleAU(double debye) { return debye * debyeToDipoleAUFactor; }
        [[nodiscard]] constexpr static double dipoleAUToDebye(double dipole) { return dipole / debyeToDipoleAUFactor; }
    };
}