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

#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"

namespace Physica {
    /**
     * \class WaterPolarTensor provides polarization tensor of water molecular at CCSD/W6 precision as introduced in [1].
     * 
     * Reference:
     * [1] J. Chem. Phys. 122, 144310 (2005); https://doi.org/10.1063/1.1867437
     */
    template<Scalar T, bool UseDynamicPolar>
    class WaterPolarTensor {
        /* Make curve pos */
        constexpr static double refBondLength = PhyConst<AU>::angstromToBohr(0.95843);
        constexpr static double refBondAngle = PhyConst<SI>::degreeToRadian(104.44);
        constexpr static double normalLength = PhyConst<AU>::angstromToBohr(1);
        /* Convert factors */
        constexpr static double UnitChargeSI = PhyConst<SI>::unitCharge;
        constexpr static double BohrSI = PhyConst<SI>::bohrRadius;
        constexpr static double HartreeSI = UnitChargeSI * PhyConst<AU>::hartreeToEv(1);
        constexpr static double Factor = UnitChargeSI * UnitChargeSI / HartreeSI;
        constexpr static double UnitDipoleToSI = Factor * BohrSI * BohrSI;
        constexpr static double convertFactor1 = 1E-40 / UnitDipoleToSI;
        constexpr static double convertFactor2 = 1E-30 / (Factor * BohrSI);
        constexpr static double convertFactor3 = 1E-20 / Factor;
        constexpr static double convertFactor4 = 1E-10 * BohrSI / Factor;
        constexpr static double convertFactor5 = BohrSI * BohrSI / Factor;

        constexpr static size_t NumDiagBase = 22;
        constexpr static size_t NumOffDiagBase = 13;
        using Matrix3D = PeriodicCell<T, 3>::LatticeMatrix;
        using DiagBaseVector = DenseVector<T, NumDiagBase>;
        using OffDiagBaseVector = DenseVector<T, NumOffDiagBase>;
    private:
        Matrix3D polarTensor; //In Molecular Frame
        Matrix3D labFrame;
    public:
        WaterPolarTensor(Vector3D<T> posOH1, Vector3D<T> posOH2);
        WaterPolarTensor(const WaterPolarTensor&) = default;
        WaterPolarTensor(WaterPolarTensor&&) noexcept = default;
        ~WaterPolarTensor() = default;
        /* Operators */
        WaterPolarTensor& operator=(WaterPolarTensor model) noexcept;
        [[nodiscard]] Vector3D<T> operator*(const Vector3D<T>& v) const;
        /* Operations */
        void swap(WaterPolarTensor& __restrict model) noexcept;
    private:
        static Matrix3D makePolarTensor(T curvePos1, T curvePos2, T curvePos3);
        static DiagBaseVector makeDiagBases(T curvePos1, T curvePos2, T curvePos3);
        static OffDiagBaseVector makeOffDiagBases(T curvePos1, T curvePos2, T curvePos3);
        static DiagBaseVector makeFactorXX();
        static DiagBaseVector makeFactorYY();
        static DiagBaseVector makeFactorZZ();
        static OffDiagBaseVector makeFactorXZ();
    };

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::WaterPolarTensor(
            Vector3D<T> posOH1, Vector3D<T> posOH2) {
        const T bondLength1 = posOH1.norm();
        const T bondLength2 = posOH2.norm();
        const T bondAngle = arccos(posOH1 * posOH2 / (bondLength1 * bondLength2));
        const T curvePos1 = (bondLength1 + bondLength2 - T(2 * refBondLength)) * T(M_SQRT1_2);
        const T curvePos2 = T(normalLength) * (bondAngle - T(refBondAngle));
        const T curvePos3 = (bondLength2 - bondLength1) * T(M_SQRT1_2);
        polarTensor = makePolarTensor(curvePos1, curvePos2, curvePos3);
        
        posOH1 *= reciprocal(bondLength1);
        posOH2 *= reciprocal(bondLength2);
        auto x = labFrame.row(0);
        x = -(posOH1 + posOH2);
        auto z = labFrame.row(2);
        z = posOH1 - posOH2;
        auto y = labFrame.row(1);
        y = cross(z, x);
        x.toUnit();
        y.toUnit();
        z.toUnit();
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>&
    WaterPolarTensor<T, UseDynamicPolar>::operator=(WaterPolarTensor model) noexcept {
        swap(model);
        return *this;
    }

    template<Scalar T, bool UseDynamicPolar>
    Vector3D<T> WaterPolarTensor<T, UseDynamicPolar>::operator*(const Vector3D<T>& v) const {
        const Vector3D<T> vInMolecularFrame = labFrame * v;
        const Vector3D<T> uInMolecularFrame = polarTensor * vInMolecularFrame;
        const Vector3D<T> u = labFrame.transpose() * uInMolecularFrame;
        return u;
    }

    template<Scalar T, bool UseDynamicPolar>
    void WaterPolarTensor<T, UseDynamicPolar>::swap(WaterPolarTensor& __restrict model) noexcept {
        assert(this != &model && "[Error]: Self swap is likely a bug");
        polarTensor.swap(model.polarTensor);
        labFrame.swap(model.labFrame);
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::Matrix3D
    WaterPolarTensor<T, UseDynamicPolar>::makePolarTensor(T curvePos1, T curvePos2, T curvePos3) {
        Matrix3D polarTensor(3, 3, 0);
        const DiagBaseVector diagBases = makeDiagBases(curvePos1, curvePos2, curvePos3);
        polarTensor[0, 0] = diagBases * makeFactorXX();
        polarTensor[1, 1] = diagBases * makeFactorYY();
        polarTensor[2, 2] = diagBases * makeFactorZZ();
        polarTensor[2, 0] = polarTensor[0, 2] = makeOffDiagBases(curvePos1, curvePos2, curvePos3) * makeFactorXZ();
        return polarTensor;
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::DiagBaseVector
    WaterPolarTensor<T, UseDynamicPolar>::makeDiagBases(T curvePos1, T curvePos2, T curvePos3) {
        DiagBaseVector bases(NumDiagBase);
        bases[0] = 1;
        bases[1] = curvePos2;
        bases[2] = curvePos1;
        bases[3] = square(curvePos3);
        bases[4] = square(curvePos2);
        bases[5] = square(curvePos1);
        bases[6] = curvePos1 + curvePos2;
        bases[7] = bases[4] * curvePos2;
        bases[8] = bases[5] * curvePos1;
        bases[9] = curvePos2 + bases[3];
        bases[10] = curvePos1 + bases[3];
        bases[11] = curvePos1 + bases[4];
        bases[12] = bases[5] + curvePos2;
        bases[13] = square(bases[3]);
        bases[14] = square(bases[4]);
        bases[15] = square(bases[5]);
        bases[16] = bases[3] + bases[4];
        bases[17] = curvePos1 + curvePos2 + bases[3];
        bases[18] = bases[3] + bases[5];
        bases[19] = bases[4] + bases[5];
        bases[20] = bases[8] + curvePos2;
        bases[21] = bases[7] + curvePos1;
        return bases;
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::OffDiagBaseVector
    WaterPolarTensor<T, UseDynamicPolar>::makeOffDiagBases(T curvePos1, T curvePos2, T curvePos3) {
        OffDiagBaseVector bases(NumOffDiagBase);
        bases[0] = curvePos3;
        bases[1] = curvePos2 + curvePos3;
        bases[2] = curvePos1 + curvePos3;
        bases[3] = square(curvePos3) * curvePos3;
        bases[4] = square(curvePos2) + curvePos3;
        bases[5] = curvePos1 + bases[1];
        bases[6] = square(curvePos1) + curvePos3;
        bases[7] = curvePos2 + bases[3];
        bases[8] = curvePos1 + bases[3];
        bases[9] = bases[2] + square(curvePos2);
        bases[10] = square(curvePos2) * curvePos2 + curvePos3;
        bases[11] = square(curvePos1) + bases[2];
        bases[12] = square(curvePos1) * curvePos1 + curvePos3;
        return bases;
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::DiagBaseVector
    WaterPolarTensor<T, UseDynamicPolar>::makeFactorXX() {
        if constexpr (UseDynamicPolar) {
            DiagBaseVector result{
                1.63133 * convertFactor1,
                -0.06044 * convertFactor2,
                1.69739 * convertFactor2,
                0.79381 * convertFactor3,
                0.14506 * convertFactor3,
                0.71438 * convertFactor3,
                -0.68859 * convertFactor3,
                0.05549 * convertFactor4,
                -0.16034 * convertFactor4,
                -1.14551 * convertFactor4,
                0.59371 * convertFactor4,
                0.02855 * convertFactor4,
                -0.31977 * convertFactor4,
                -0.18252 * convertFactor5,
                -0.03486 * convertFactor5,
                -0.54566 * convertFactor5,
                0.37602 * convertFactor5,
                -0.03442 * convertFactor5,
                -1.02788 * convertFactor5,
                0.01697 * convertFactor5,
                0.24590 * convertFactor5,
                0.18296 * convertFactor5};
            return result;
        }
        else {
            DiagBaseVector result{
                1.58787 * convertFactor1,
                -0.07614 * convertFactor2,
                1.60710 * convertFactor2,
                0.71119 * convertFactor3,
                0.13331 * convertFactor3,
                0.61842 * convertFactor3,
                -0.70521 * convertFactor3,
                0.05328 * convertFactor4,
                -0.20959 * convertFactor4,
                -1.08201 * convertFactor4,
                0.35624 * convertFactor4,
                0.00892 * convertFactor4,
                -0.35495 * convertFactor4,
                -0.22179 * convertFactor5,
                -0.03182 * convertFactor5,
                -0.51709 * convertFactor5,
                0.34644 * convertFactor5,
                0.01374 * convertFactor5,
                -1.25979 * convertFactor5,
                -0.03825 * convertFactor5,
                0.18326 * convertFactor5,
                0.17426 * convertFactor5};
            return result;
        }
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::DiagBaseVector
    WaterPolarTensor<T, UseDynamicPolar>::makeFactorYY() {
        if constexpr (UseDynamicPolar) {
            DiagBaseVector result{
                1.59884 * convertFactor1,
                0.08866 * convertFactor2,
                0.90706 * convertFactor2,
                -0.03129 * convertFactor3,
                0.06210 * convertFactor3,
                0.14658 * convertFactor3,
                0.10276 * convertFactor3,
                -0.00325 * convertFactor4,
                -0.06520 * convertFactor4,
                -0.10133 * convertFactor4,
                -0.04683 * convertFactor4,
                0.10669 * convertFactor4,
                0.11117 * convertFactor4,
                -0.05699 * convertFactor5,
                -0.00847 * convertFactor5,
                -0.06240 * convertFactor5,
                0.17064 * convertFactor5,
                -0.17870 * convertFactor5,
                -0.34393 * convertFactor5,
                0.02100 * convertFactor5,
                -0.02978 * convertFactor5,
                0.06393 * convertFactor5};
            return result;
        }
        else {
            DiagBaseVector result{
                1.53818 * convertFactor1,
                0.08361 * convertFactor2,
                0.82863 * convertFactor2,
                -0.04232 * convertFactor3,
                0.05666 * convertFactor3,
                0.08172 * convertFactor3,
                0.08557 * convertFactor3,
                -0.00271 * convertFactor4,
                -0.11892 * convertFactor4,
                -0.06926 * convertFactor4,
                -0.10806 * convertFactor4,
                0.08190 * convertFactor4,
                0.06999 * convertFactor4,
                -0.05890 * convertFactor5,
                -0.00800 * convertFactor5,
                -0.11111 * convertFactor5,
                0.12995 * convertFactor5,
                -0.09119 * convertFactor5,
                -0.48361 * convertFactor5,
                -0.03088 * convertFactor5,
                -0.09411 * convertFactor5,
                0.05992 * convertFactor5};
            return result;
        }
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::DiagBaseVector
    WaterPolarTensor<T, UseDynamicPolar>::makeFactorZZ() {
        if constexpr (UseDynamicPolar) {
            DiagBaseVector result{
                1.68167 * convertFactor1,
                0.25330 * convertFactor2,
                2.49618 * convertFactor2,
                1.00569 * convertFactor3,
                0.04642 * convertFactor3,
                2.17892 * convertFactor3,
                1.03261 * convertFactor3,
                -0.09462 * convertFactor4,
                0.73744 * convertFactor4,
                1.09418 * convertFactor4,
                0.76823 * convertFactor4,
                -0.06618 * convertFactor4,
                1.73783 * convertFactor4,
                0.11518 * convertFactor5,
                -0.02222 * convertFactor5,
                -0.51425 * convertFactor5,
                -0.16139 * convertFactor5,
                4.05248 * convertFactor5,
                -1.89819 * convertFactor5,
                0.21050 * convertFactor5,
                1.41546 * convertFactor5,
                -0.20008 * convertFactor5};
            return result;
        }
        else {
            DiagBaseVector result{
                1.64720 * convertFactor1,
                0.24730 * convertFactor2,
                2.39394 * convertFactor2,
                0.89549 * convertFactor3,
                0.03954 * convertFactor3,
                1.99483 * convertFactor3,
                0.99680 * convertFactor3,
                -0.08900 * convertFactor4,
                0.53619 * convertFactor4,
                0.96599 * convertFactor4,
                0.41179 * convertFactor4,
                -0.08290 * convertFactor4,
                1.62467 * convertFactor4,
                0.04290 * convertFactor5,
                -0.02203 * convertFactor5,
                -0.62916 * convertFactor5,
                -0.20222 * convertFactor5,
                3.29320 * convertFactor5,
                -2.27559 * convertFactor5,
                0.15971 * convertFactor5,
                1.20526 * convertFactor5,
                -0.18405 * convertFactor5};
            return result;
        }
    }

    template<Scalar T, bool UseDynamicPolar>
    WaterPolarTensor<T, UseDynamicPolar>::OffDiagBaseVector
    WaterPolarTensor<T, UseDynamicPolar>::makeFactorXZ() {
        if constexpr (UseDynamicPolar) {
            OffDiagBaseVector result{
                1.24322 * convertFactor2,
                -0.36918 * convertFactor3,
                2.06801 * convertFactor3,
                0.08804 * convertFactor4,
                -0.35436 * convertFactor4,
                0.08085 * convertFactor4,
                0.93518 * convertFactor4,
                0.60596 * convertFactor5,
                -0.44015 * convertFactor5,
                -0.58108 * convertFactor5,
                -0.14033 * convertFactor5,
                1.71445 * convertFactor5,
                -1.08751 * convertFactor5};
            return result;
        }
        else {
            OffDiagBaseVector result{
                1.17904 * convertFactor2,
                -0.37700 * convertFactor3,
                1.85559 * convertFactor3,
                0.00698 * convertFactor4,
                -0.36099 * convertFactor4,
                -0.00202 * convertFactor4,
                0.62070 * convertFactor4,
                0.53070 * convertFactor5,
                -0.64252 * convertFactor5,
                -0.63496 * convertFactor5,
                -0.12823 * convertFactor5,
                1.40356 * convertFactor5,
                -1.24044 * convertFactor5};
            return result;
        }
    }
}
