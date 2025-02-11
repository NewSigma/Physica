/*
 * Copyright 2021-2025 Weibo He.
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
#include <numbers>

namespace Physica {
    enum UnitSystem {
        SI,
        AU,
        ESU,
        QE
    };

    template<UnitSystem unitSystem = SI> class PhyConst;
    /**
     * Physical constants are from [1], which updates every 4 years.
     * Relative atom mass are from [2]
     * 
     * Note:
     * Relative atom mass are taken to be the average value of isotopes that are common in nature as [2] does. If you are interesed in
     * a particular isotope, modifications are necessary.
     * 
     * Reference:
     * [1] Physical constants, SI (NIST 2022); http://physics.nist.gov/constants
     * [2] Standard Atomic Weights 2021; https://www.ciaaw.org/atomic-weights.htm
     */
    template<>
    class PhyConst<SI> {
        constexpr static double calorieInJoule = 4.184; // Reference: http://simple.wikipedia.org/wiki/Calorie
    public:
        /* Universal */
        constexpr static double gravitationG = 6.67430E-11;
        constexpr static double planck = 6.62607015E-34;
        constexpr static double reducedPlanck = planck / (2 * std::numbers::pi);
        constexpr static double speedOfLight = 299792458;
        constexpr static double vacuumDielectric = 8.8541878188E-12;
        /* Electromagnetic */
        constexpr static double unitCharge = 1.602176634E-19;
        /* Atomic and nuclear */
        constexpr static double bohrRadius = 0.529177210544E-10;
        constexpr static double electronMass = 9.1093837139E-31;
        constexpr static double protonMass = 1.67262192595E-27;
        constexpr static double neutronMass = 1.67492750056E-27;
        /* Physico-chemical */
        constexpr static double atomicMassConst = 1.66053906892E-27;
        constexpr static double avogadroNa = 6.02214076E23;
        constexpr static double boltzmannK = 1.380649E-23;
        /**
         * The first element is a space holder
         */
        constexpr static double relativeAtomMass[37]{0, 1.00798, 4.002602,
                6.968, 9.0121831, 10.814, 12.0106, 14.00686, 15.99940, 18.998403162, 20.1797,
                22.98976928, 24.306, 26.9815384, 28.085, 30.973761998, 32.068, 35.452, 39.878,
                39.0983, 40.078, 44.955907, 47.867, 50.9415, 51.9961, 54.938043, 55.845, 58.933194, 58.6934, 63.546, 65.38, 69.723, 72.630, 74.921595, 78.971, 79.904, 83.798};
        constexpr static const char* elementSymbol[37]{"-", "H", "He",
                "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"};
    
        [[nodiscard]] constexpr static double degreeToRadian(double degree) { return degree / 180 * M_PI; }
        [[nodiscard]] constexpr static double radianToDegree(double radian) { return radian / M_PI * 180; }
        [[nodiscard]] constexpr static double atomMass(size_t atomicNum) { return relativeAtomMass[atomicNum] * atomicMassConst; }
        [[nodiscard]] constexpr static double calorieToJoule(double calorie) { return calorie * calorieInJoule; }
        [[nodiscard]] constexpr static double waveNumToTHz(double waveNum) { return waveNum * speedOfLight * 1E-10; }
    };
    /**
     * Hartree atomic units. \param hartreeInEv referenced from [1]
     * 
     * Reference:
     * [1] Physical constants, SI (NIST 2022); http://physics.nist.gov/constants
     */
    template<>
    class PhyConst<AU> {
    private:
        constexpr static double hartreeInEv = 27.211386245981;
        constexpr static double rydbergInEv = hartreeInEv * 0.5;
        constexpr static double bohrInAngstorm = PhyConst<SI>::bohrRadius * 1E10;
        constexpr static double hartreeInJoule = PhyConst<SI>::unitCharge * hartreeInEv;
        constexpr static double jouleInHartree = 1 / hartreeInJoule;
        constexpr static double timeInSecond = PhyConst<SI>::reducedPlanck * jouleInHartree;
        constexpr static double temperatureInK = 1 / (PhyConst<SI>::boltzmannK * jouleInHartree);
        constexpr static double pressInGPa = hartreeInJoule / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius) * 1E-9;
    public:
        /* Universal */
        constexpr static double planck = M_PI * 2;
        constexpr static double reducedPlanck = 1;
        constexpr static double speedOfLight = PhyConst<SI>::speedOfLight / PhyConst<SI>::bohrRadius * timeInSecond;
        constexpr static double vacuumDielectric = (PhyConst<SI>::vacuumDielectric * PhyConst<SI>::bohrRadius) / (jouleInHartree * PhyConst<SI>::unitCharge * PhyConst<SI>::unitCharge);
        /* Electromagnetic */
        constexpr static double unitCharge = 1;
        /* Atomic and nuclear */
        constexpr static double bohrRadius = 1;
        constexpr static double electronMass = 1;
        constexpr static double protonMass = PhyConst<SI>::protonMass / PhyConst<SI>::electronMass;
        constexpr static double neutronMass = PhyConst<SI>::neutronMass / PhyConst<SI>::electronMass;
        /* Physico-chemical */
        constexpr static double boltzmannK = 1;

        [[nodiscard]] constexpr static double hartreeToEv(double hartree) { return hartree * hartreeInEv; }
        [[nodiscard]] constexpr static double eVToHartree(double ev) { return ev * (1.0 / hartreeInEv); }
        [[nodiscard]] constexpr static double bohrToAngstorm(double bohr) { return bohr * bohrInAngstorm; }
        [[nodiscard]] constexpr static double angstormToBohr(double angstorm) { return angstorm * (1.0 / bohrInAngstorm); }
        [[nodiscard]] constexpr static double timeToSecond(double atomic_time) { return atomic_time * timeInSecond; }
        [[nodiscard]] constexpr static double secondToTime(double second) { return second / timeInSecond; }
        [[nodiscard]] constexpr static double temperatureToK(double atomic_tem) { return atomic_tem * temperatureInK; }
        [[nodiscard]] constexpr static double kToTemperature(double kelvin) { return kelvin / temperatureInK; }
        [[nodiscard]] constexpr static double pressToGPa(double atomic_press) { return atomic_press * pressInGPa; }
        [[nodiscard]] constexpr static double pressToKBar(double atomic_press) { return pressToGPa(atomic_press) * 10; }
        [[nodiscard]] constexpr static double atomMass(size_t atomicNum) { return PhyConst<SI>::atomMass(atomicNum) / PhyConst<SI>::electronMass; }
        [[nodiscard]] constexpr static double freqToTHz(double atomic_freq) { return atomic_freq * 1E-12 / timeToSecond(1); }
        [[nodiscard]] constexpr static double THzToFreq(double thz) { return thz * timeToSecond(1) * 1E12; }
    };

    template<>
    class PhyConst<ESU> {
        constexpr static double esuToCoulombFactor = 3.335640951084828E-10; // std::sqrt(PhyConst<SI>::vacuumDielectric * 4 * M_PI * 10E-9);
    public:
        constexpr static double unitCharge = PhyConst<SI>::unitCharge / esuToCoulombFactor;
    private:
        constexpr static double debyeToDipoleAUFactor = (1E-10 / unitCharge) * PhyConst<AU>::angstormToBohr(1);
    public:
        [[nodiscard]] constexpr static double esuToCoulomb(double esu) { return esu * esuToCoulombFactor; }
        [[nodiscard]] constexpr static double debyeToDipoleAU(double debye) { return debye * debyeToDipoleAUFactor; }
        [[nodiscard]] constexpr static double dipoleAUToDebye(double dipole) { return dipole / debyeToDipoleAUFactor; }
    };

    template<>
    class PhyConst<QE> {
    public:
        constexpr static double planck = M_PI * 4;
        constexpr static double reducedPlanck = 2;
    };
}