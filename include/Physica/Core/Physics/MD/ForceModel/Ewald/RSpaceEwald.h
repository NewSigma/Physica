/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"
#include "Physica/Core/Parallel/Parallel.h"

namespace Physica {
    /**
     * Reference:
     * [1] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications[M]; San Diego: Academic, 2002:304-306
     */
    template<Scalar T, bool IsSmallCell = false>
    class RSpaceEwald : public PairModel<RSpaceEwald<T, IsSmallCell>> {
        using This = RSpaceEwald<T, IsSmallCell>;
        using Base = PairModel<This>;
    public:
        using ComplexType = Complex<T>;
        using Base::Dim;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Vec3D;
        using typename Base::ForceConstMatrix;
        using SearchRangeType = PeriodicCell<T, Dim>::SearchRangeType;
        using Matrix3D = DenseMatrix<T, MatrixMajor::Col, 3, 3>;
        using BornChargeArray = Array<Matrix3D>;
        constexpr static size_t ErfcTableSize = 4096 + 512 + 2;
        constexpr static double ErfcTableStep = 0.001;
        constexpr static double SumPrec = (ErfcTableSize - 2) * ErfcTableStep; // Referenced from [1], minus 2 to avoid overflow
    protected:
        using typename Base::Tv;
    private:
        LatticeMatrix lattice;
        LatticeMatrix repLatt;
        VectorND<T> charges;
        VectorND<T> erfc_table;
        CoDiff<T> volume;
        CoDiff<T> inv_volume;
        Tv integralLimit;
        Tv erfcStep;
        Tv repErfcStep;
        Tv repDoubleSquareStep;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
    public:
        RSpaceEwald() = default;
        RSpaceEwald(LatticeMatrix lattice_, VectorND<T> charges_);
        RSpaceEwald(const RSpaceEwald&) = default;
        RSpaceEwald(RSpaceEwald&&) noexcept = default;
        ~RSpaceEwald() = default;
        /* Operator */
        RSpaceEwald& operator=(RSpaceEwald obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T potentialV(const PositionMatrix& pos) const;

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const PositionMatrix& pos) const;

        [[nodiscard]] ComplexType forceConst(const PositionMatrix& pos, const Vec3D& waveQ, size_t dof1, size_t dof2) const;

        [[nodiscard]] LatticeMatrix virial(const PositionMatrix& pos) const;
        [[nodiscard]] BornChargeArray calcBornCharge() const { return makeBornCharge(charges); }
        void swap(RSpaceEwald& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const auto& getRepLattice() const noexcept { return repLatt; }
        [[nodiscard]] const auto& getCharges() const noexcept { return charges; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
        [[nodiscard]] T getVolume() const noexcept { return volume; }
        [[nodiscard]] T getInvVolume() const noexcept { return inv_volume; }
        [[nodiscard]] Tv getIntegralLimit() const noexcept { return integralLimit; }
        [[nodiscard]] T getRSpaceCutoff() const noexcept { return Base::getCutoff(); }
        [[nodiscard]] T getSquaredRSpaceCutoff() const noexcept { return Base::getSquaredCutoff(); }
        [[nodiscard]] const auto& getRSpaceSumRange() const noexcept { return rSpaceSumRange; }
        [[nodiscard]] const auto& getKSpaceSumRange() const noexcept { return kSpaceSumRange; }
        /* Setters */
        void setLattice(LatticeMatrix lattice_);
        void setIntegralLimit(Tv integralLimit_);
    protected:
        [[nodiscard]] CoDiff<T> calcSelfE() const;
        [[nodiscard]] CoDiff<T> calcGammaPointE() const;
        [[nodiscard]] CoDiff<T> pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] CoDiff<T> force_functor(size_t i, size_t j, T r, T r2) const;
    private:
        /* Operations */
        void makeTables();
        [[nodiscard]] T pot_functor_slow(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] T force_functor_slow(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] T rSpaceForceConstImpl1(T r) const;
        [[nodiscard]] T rSpaceForceConstImpl2(T r) const;
        using Base::potentialV;
        using Base::forceConst;
        /* Getters */
        using Base::getCutoff;
        using Base::getSquaredCutoff;
        /* Static members */
        [[nodiscard]] static BornChargeArray makeBornCharge(const VectorND<T>& charges);
        /* Friends */
        friend class PairModel<RSpaceEwald<T, IsSmallCell>>;
        friend class device_obj<This>;
        friend class Physica::Test;
    };

    template<Scalar T, bool IsSmallCell>
    RSpaceEwald<T, IsSmallCell>::RSpaceEwald(LatticeMatrix lattice_, VectorND<T> charges_)
            : charges(std::move(charges_)), erfc_table(ErfcTableSize + 1) {
        setLattice(std::move(lattice_));
    }
    /**
     * \param pos must be in cartesian convention
     */
    template<Scalar T, bool IsSmallCell>
    T RSpaceEwald<T, IsSmallCell>::potentialV(const PositionMatrix& pos) const {
        return Base::potentialV(lattice, pos);
    }

    template<Scalar T, bool IsSmallCell>
    template<ExecutePolicy P> 
    VectorND<T> RSpaceEwald<T, IsSmallCell>::force_short(const PositionMatrix& pos) const {
        const VectorND<T> rSpaceSum = Base::template force<Sequential>(lattice, pos);
        return rSpaceSum;
    }
    /**
     * Reference:
     * [1] Rev. Mod. Phys. 73, 515; https://doi.org/10.1103/RevModPhys.73.515
     */
    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::forceConst(const PositionMatrix& pos, const Vec3D& waveQ, size_t dof1, size_t dof2) const -> ComplexType {
        const size_t atom1 = dof1 / 3;
        const size_t atom2 = dof2 / 3;
        const size_t direction1 = dof1 % 3U;
        const size_t direction2 = dof2 % 3U;
        const T charge1 = charges[atom1];
        const T charge2 = charges[atom2];

        ComplexType rSpaceSum = 0;
        PeriodicCell<T, Dim>::forCellInRange(rSpaceSumRange, lattice,
            [this, atom1, atom2, direction1, direction2, charge1, charge2, waveQ, &pos, &rSpaceSum](Vec3D delta) {
                const bool isSameDirection = direction1 == direction2;
                ComplexType temp = 0;
                {
                    const Vec3D x = pos.row(atom1) - pos.row(atom2) + delta;
                    const T norm = x.norm();
                    const bool isNotSelf = !norm.isSubNormal();
                    if (isNotSelf) {
                        const T phase = waveQ * delta;
                        temp = rSpaceForceConstImpl1(norm) * (x[direction1] * x[direction2]);
                        if (isSameDirection)
                            temp += rSpaceForceConstImpl2(norm);
                        temp *= ComplexType::fromPhase(phase) * (charge1 * charge2);
                    }
                }

                const bool isSameAtom = atom1 == atom2;
                if (isSameAtom) {
                    ComplexType temp1 = 0;
                    for (size_t i = 0; i < getNumParticle(); ++i) {
                        const Vec3D x = pos.row(atom1) - pos.row(i) + delta;
                        const T norm = x.norm();
                        const bool isNotSelf = !norm.isSubNormal();
                        if (isNotSelf) {
                            ComplexType temp2 = rSpaceForceConstImpl1(norm) * (x[direction1] * x[direction2]);
                            if (isSameDirection)
                                temp2 += rSpaceForceConstImpl2(norm);
                            temp1 += temp2 * charges[i];
                        }
                    }
                    temp -= temp1 * charge1;
                }
                rSpaceSum += temp;
            });
        return rSpaceSum;
    }

    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::virial(const PositionMatrix& pos) const -> LatticeMatrix {
        return Base::virial(lattice, pos);
    }

    template<Scalar T, bool IsSmallCell>
    void RSpaceEwald<T, IsSmallCell>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        lattice.swap(obj.lattice);
        repLatt.swap(obj.repLatt);
        charges.swap(obj.charges);
        erfc_table.swap(obj.erfc_table);
        volume.swap(obj.volume);
        inv_volume.swap(obj.inv_volume);
        integralLimit.swap(obj.integralLimit);
        erfcStep.swap(obj.erfcStep);
        repErfcStep.swap(obj.repErfcStep);
        repDoubleSquareStep.swap(obj.repDoubleSquareStep);
        rSpaceSumRange.swap(obj.rSpaceSumRange);
        kSpaceSumRange.swap(obj.kSpaceSumRange);
    }

    template<Scalar T, bool IsSmallCell>
    void RSpaceEwald<T, IsSmallCell>::setLattice(LatticeMatrix lattice_) {
        assert(charges.getLength() != 0 && "[Error]: Charges should be initialized before lattice update");
        lattice = std::move(lattice_);
        repLatt = PeriodicCell<T, Dim>::makeRepLattice(lattice);
        volume = PeriodicCell<T, Dim>::getVolume(lattice);
        inv_volume = reciprocal(volume);

        const Tv averageCellSize = cbrt(volume.value());
        const Tv estimate = sqrt(cbrt(Tv(getNumParticle())) * Tv(M_PI)) / averageCellSize;
        setIntegralLimit(estimate);
    }

    template<Scalar T, bool IsSmallCell>
    void RSpaceEwald<T, IsSmallCell>::setIntegralLimit(Tv integralLimit_) {
        assert(integralLimit_.isPositive() && "[Error]: Invalid integralLimit");
        const Tv heightX_2Pi = reciprocal(repLatt.row(0).values().norm());
        const Tv heightY_2Pi = reciprocal(repLatt.row(1).values().norm());
        const Tv heightZ_2Pi = reciprocal(repLatt.row(2).values().norm());
        const Tv factor1 = 2 * M_PI * (1 - std::numeric_limits<T>::epsilon()); //To avoid rSpaceCutoff larger than max value
        const Tv maxRSpaceCutoff = std::min(heightX_2Pi, std::min(heightY_2Pi, heightZ_2Pi)) * Tv(factor1);
        const Tv minLimit = Tv(SumPrec) / maxRSpaceCutoff;
        integralLimit = std::max(integralLimit_, minLimit);

        const Tv rSpaceCutoff = Tv(SumPrec) / integralLimit;
        rSpaceSumRange = PeriodicCell<T, Dim>::estimateRange(lattice, rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<T, Dim>::estimateRange(repLatt, Tv(SumPrec * 2) * integralLimit);
        makeTables();
        Base::setCutoff(rSpaceCutoff);
    }

    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::calcSelfE() const -> CoDiff<T> {
        return square(charges).sum() * (integralLimit / sqrt(Tv(M_PI)));
    }

    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::calcGammaPointE() const -> CoDiff<T> {
        return square(charges.sum()) * Tv(-M_PI) / (Tv(2) * square(integralLimit)) * inv_volume;
    }
    /**
     * Optimize: make use of x1, x2, x3 are equal distance
     */
    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::pot_functor(size_t i, size_t j, T r, T) const -> CoDiff<T> {
        const Tv temp = r.value() * repErfcStep + Tv(0.5);
        const int index = double(temp);
        const Tv x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        const T interp = Internal::quadraticInterpolate<T>(x1 - erfcStep, x1, x1 + erfcStep, y[0], y[1], y[2], r);
        return charges[i] * charges[j] * interp;
    }

    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::force_functor(size_t i, size_t j, T r, T) const -> CoDiff<T> {
        const Tv temp = r.value() * repErfcStep + Tv(0.5);
        const int index = double(temp);
        const Tv x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        return -charges[i] * charges[j] * Internal::quadraticInterpolate_diff1<T>(repDoubleSquareStep, erfcStep, x1, y[0], y[1], y[2], r);
    }

    template<Scalar T, bool IsSmallCell>
    void RSpaceEwald<T, IsSmallCell>::makeTables() {
        for (size_t i = 2; i < erfc_table.getLength(); ++i) {
            const auto x = Tv(ErfcTableStep * double(i - 1));
            erfc_table[i] = erfc(x) / x * integralLimit;
        }
        erfc_table[0] = erfc_table[1] = erfc_table[2]; // Smooth out divergent erfc(0) / 0 
        erfcStep = Tv(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        repDoubleSquareStep = reciprocal(square(erfcStep) * Tv(2));
    }
    /**
     * Slow version functors are provided for debug use
     */
    template<Scalar T, bool IsSmallCell>
    T RSpaceEwald<T, IsSmallCell>::pot_functor_slow(size_t i, size_t j, T r, [[maybe_unused]] T r2) const {
        const T x = r * integralLimit;
        return charges[i] * charges[j] * erfc(x) / r;
    }

    template<Scalar T, bool IsSmallCell>
    T RSpaceEwald<T, IsSmallCell>::force_functor_slow(size_t i, size_t j, T r, T r2) const {
        const T x = r * integralLimit;
        return charges[i] * charges[j] * (erfc(x) + x * exp(-square(x)) * Tv(M_2_SQRTPI)) / r2;
    }

    template<Scalar T, bool IsSmallCell>
    T RSpaceEwald<T, IsSmallCell>::rSpaceForceConstImpl1(T r) const {
        const T x = integralLimit * r;
        const T x2 = square(x);
        const T term1 = T(3) * erfc(x);
        const T term2 = T(2 * M_2_SQRTPI) * x * (T(1.5) + x2) / exp(x2);
        return (term1 + term2) / (square(square(r)) * r);
    }

    template<Scalar T, bool IsSmallCell>
    T RSpaceEwald<T, IsSmallCell>::rSpaceForceConstImpl2(T r) const {
        const T x = integralLimit * r;
        const T term1 = erfc(x);
        const T term2 = T(M_2_SQRTPI) * x / exp(square(x));
        return -(term1 + term2) / (square(r) * r);
    }

    template<Scalar T, bool IsSmallCell>
    auto RSpaceEwald<T, IsSmallCell>::makeBornCharge(const VectorND<T>& charges) -> BornChargeArray {
        BornChargeArray result(charges.getLength(), Matrix3D(3, 3, T(0)));
        for (size_t i = 0; i < result.getLength(); ++i) {
            auto diag = result[i].diag();
            diag = charges[i];
        }
        return result;
    }
}

namespace Physica {
    template<class T, bool B>
    class Traits<RSpaceEwald<T, B>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPotDependOnAtomIndex = true;
        constexpr static bool IsSmallCell = B;
    };
}
