/*
 * Copyright 2021-2024 Weibo He.
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
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] D. Frenkel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications; San Diego: Academic, 2002:304-306
     */
    template<class ScalarType, bool IsSmallCell = false>
    class RSpaceEwald : public PairModel<RSpaceEwald<ScalarType, IsSmallCell>> {
        using This = RSpaceEwald<ScalarType, IsSmallCell>;
        using Base = PairModel<This>;
    public:
        using ComplexType = Complex<ScalarType>;
        using Base::Dim;
        using typename Base::PlainScalar;
        using typename Base::LatticeMatrix;
        using typename Base::PositionMatrix;
        using typename Base::Vector3D;
        using typename Base::ForceConstMatrix;
        using SearchRangeType = typename PeriodicCell<ScalarType, Dim>::SearchRangeType;
        using Matrix3D = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector, 3, 3>;
        using BornChargeArray = Array<Matrix3D>;
        constexpr static size_t ErfcTableSize = 4096 + 512 + 2;
        constexpr static double ErfcTableStep = 0.001;
        constexpr static double SumPrec = (ErfcTableSize - 2) * ErfcTableStep; // Referenced from [1], minus 2 to avoid overflow
    private:
        LatticeMatrix lattice;
        LatticeMatrix repLatt;
        Vector<ScalarType> charges;
        Vector<ScalarType> erfc_table;
        ScalarType volume;
        ScalarType inv_volume;
        ScalarType integralLimit;
        ScalarType erfcStep;
        ScalarType repErfcStep;
        ScalarType repDoubleSquareStep;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
    public:
        RSpaceEwald() = default;
        RSpaceEwald(LatticeMatrix lattice_, Vector<ScalarType> charges_);
        RSpaceEwald(const RSpaceEwald&) = default;
        RSpaceEwald(RSpaceEwald&&) noexcept = default;
        ~RSpaceEwald() = default;
        /* Operator */
        RSpaceEwald& operator=(RSpaceEwald obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] inline ScalarType potentialV(const PositionMatrix& pos) const;

        template<class Executor>
        [[nodiscard]] inline Vector<ScalarType> force_short(const PositionMatrix& pos) const;

        [[nodiscard]] ComplexType forceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const;

        [[nodiscard]] inline LatticeMatrix virial(const PositionMatrix& pos) const;
        [[nodiscard]] BornChargeArray calcBornCharge() const { return makeBornCharge(charges); }
        void swap(RSpaceEwald& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const LatticeMatrix& getRepLattice() const noexcept { return repLatt; }
        [[nodiscard]] const Vector<ScalarType>& getCharges() const noexcept { return charges; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return charges.getLength(); }
        [[nodiscard]] ScalarType getVolume() const noexcept { return volume; }
        [[nodiscard]] ScalarType getInvVolume() const noexcept { return inv_volume; }
        [[nodiscard]] ScalarType getIntegralLimit() const noexcept { return integralLimit; }
        [[nodiscard]] ScalarType getRSpaceCutoff() const noexcept { return Base::getCutoff(); }
        [[nodiscard]] ScalarType getSquaredRSpaceCutoff() const noexcept { return Base::getSquaredCutoff(); }
        [[nodiscard]] const SearchRangeType& getRSpaceSumRange() const noexcept { return rSpaceSumRange; }
        [[nodiscard]] const SearchRangeType& getKSpaceSumRange() const noexcept { return kSpaceSumRange; }
        /* Setters */
        void setLattice(LatticeMatrix lattice_);
        void setIntegralLimit(ScalarType integralLimit_);
    protected:
        [[nodiscard]] inline ScalarType calcSelfE() const;
        [[nodiscard]] inline ScalarType calcGammaPointE() const;
        [[nodiscard]] inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
    private:
        /* Operations */
        void makeTables();
        [[nodiscard]] inline ScalarType pot_functor_slow(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] inline ScalarType force_functor_slow(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] ScalarType rSpaceForceConstImpl1(ScalarType r) const;
        [[nodiscard]] ScalarType rSpaceForceConstImpl2(ScalarType r) const;
        using Base::potentialV;
        using Base::forceConst;
        /* Getters */
        using Base::getCutoff;
        using Base::getSquaredCutoff;
        /* Static members */
        [[nodiscard]] static BornChargeArray makeBornCharge(const Vector<ScalarType>& charges);
        /* Friends */
        friend class PairModel<RSpaceEwald<ScalarType, IsSmallCell>>;
        friend class device_obj<This>;
        friend class Physica::Test;
    };

    template<class ScalarType, bool IsSmallCell>
    RSpaceEwald<ScalarType, IsSmallCell>::RSpaceEwald(LatticeMatrix lattice_, Vector<ScalarType> charges_)
            : charges(std::move(charges_)), erfc_table(ErfcTableSize + 1) {
        setLattice(std::move(lattice_));
    }
    /**
     * \param pos must be in cartesian convension
     */
    template<class ScalarType, bool IsSmallCell>
    inline ScalarType RSpaceEwald<ScalarType, IsSmallCell>::potentialV(const PositionMatrix& pos) const {
        return Base::potentialV(lattice, pos);
    }

    template<class ScalarType, bool IsSmallCell>
    template<class Executor> 
    inline Vector<ScalarType> RSpaceEwald<ScalarType, IsSmallCell>::force_short(const PositionMatrix& pos) const {
        const Vector<ScalarType> rSpaceSum = Base::template force<SequentialExecutor>(lattice, pos);
        return rSpaceSum;
    }
    /**
     * Reference:
     * [1] Rev. Mod. Phys. 73, 515; https://doi.org/10.1103/RevModPhys.73.515
     */
    template<class ScalarType, bool IsSmallCell>
    typename RSpaceEwald<ScalarType, IsSmallCell>::ComplexType
    RSpaceEwald<ScalarType, IsSmallCell>::forceConst(const PositionMatrix& pos, const Vector3D& waveQ, size_t dof1, size_t dof2) const {
        const size_t atom1 = dof1 / 3;
        const size_t atom2 = dof2 / 3;
        const size_t direction1 = dof1 % 3U;
        const size_t direction2 = dof2 % 3U;
        const ScalarType charge1 = charges[atom1];
        const ScalarType charge2 = charges[atom2];

        ComplexType rSpaceSum = 0;
        PeriodicCell<ScalarType, Dim>::forCellInRange(rSpaceSumRange, lattice,
            [this, atom1, atom2, direction1, direction2, charge1, charge2, waveQ, &pos, &rSpaceSum](Vector3D delta) {
                const bool isSameDirection = direction1 == direction2;
                ComplexType temp = 0;
                {
                    const Vector3D x = pos.row(atom1).asVector() - pos.row(atom2).asVector() + delta;
                    const ScalarType norm = x.norm();
                    const bool isNotSelf = norm > ScalarType(std::numeric_limits<ScalarType>::min());
                    if (isNotSelf) {
                        const ScalarType phase = waveQ * delta;
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
                        const Vector3D x = pos.row(atom1).asVector() - pos.row(i).asVector() + delta;
                        const ScalarType norm = x.norm();
                        const bool isNotSelf = norm > ScalarType(std::numeric_limits<ScalarType>::min());
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

    template<class ScalarType, bool IsSmallCell>
    inline typename RSpaceEwald<ScalarType, IsSmallCell>::LatticeMatrix
    RSpaceEwald<ScalarType, IsSmallCell>::virial(const PositionMatrix& pos) const {
        return Base::virial(lattice, pos);
    }

    template<class ScalarType, bool IsSmallCell>
    void RSpaceEwald<ScalarType, IsSmallCell>::swap(RSpaceEwald& __restrict obj) noexcept {
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

    template<class ScalarType, bool IsSmallCell>
    void RSpaceEwald<ScalarType, IsSmallCell>::setLattice(LatticeMatrix lattice_) {
        assert(charges.getLength() != 0 && "[Error]: Charges should be initialized before lattice update");
        lattice = std::move(lattice_);
        repLatt = PeriodicCell<ScalarType, Dim>::makeRepLattice(lattice);
        volume = PeriodicCell<ScalarType, Dim>::getVolume(lattice);
        inv_volume = reciprocal(volume);

        const ScalarType averageCellSize = cbrt(ScalarType(volume));
        const ScalarType estimate = sqrt(cbrt(ScalarType(getNumParticle())) * ScalarType(M_PI)) / averageCellSize;
        setIntegralLimit(estimate);
    }

    template<class ScalarType, bool IsSmallCell>
    void RSpaceEwald<ScalarType, IsSmallCell>::setIntegralLimit(ScalarType integralLimit_) {
        assert(integralLimit_.isPositive() && "[Error]: Invalid integralLimit");
        const ScalarType heightX_2Pi = reciprocal(repLatt.row(0).norm());
        const ScalarType heightY_2Pi = reciprocal(repLatt.row(1).norm());
        const ScalarType heightZ_2Pi = reciprocal(repLatt.row(2).norm());
        constexpr double factor1 = 2 * M_PI * (1 - std::numeric_limits<ScalarType>::epsilon()); //To avoid rSpaceCutoff larger than max value
        const ScalarType maxRSpaceCutoff = std::min(heightX_2Pi, std::min(heightY_2Pi, heightZ_2Pi)) * PlainScalar(factor1);
        const ScalarType minLimit = ScalarType(SumPrec) / maxRSpaceCutoff;
        integralLimit = std::max(integralLimit_, minLimit).getValue();

        const PlainScalar rSpaceCutoff = PlainScalar(SumPrec) / integralLimit.getValue();
        rSpaceSumRange = PeriodicCell<ScalarType, Dim>::estimateRange(lattice, rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<ScalarType, Dim>::estimateRange(repLatt, PlainScalar(SumPrec * 2) * integralLimit.getValue());
        makeTables();
        Base::setCutoff(rSpaceCutoff);
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType RSpaceEwald<ScalarType, IsSmallCell>::calcSelfE() const {
        return square(charges).sum() * integralLimit / sqrt(PlainScalar(M_PI));
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType RSpaceEwald<ScalarType, IsSmallCell>::calcGammaPointE() const {
        return square(charges.sum()) * PlainScalar(-M_PI) / (PlainScalar(2) * square(integralLimit)) * inv_volume;
    }
    /**
     * Optimize: make use of x1, x2, x3 are equal distance
     */
    template<class ScalarType, bool IsSmallCell>
    inline ScalarType RSpaceEwald<ScalarType, IsSmallCell>::pot_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = double(temp);
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        const ScalarType interp = Internal::quadraticInterpolate<ScalarType>(x1 - erfcStep, x1, x1 + erfcStep, y[0], y[1], y[2], r);
        return charges[i] * charges[j] * interp;
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType RSpaceEwald<ScalarType, IsSmallCell>::force_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = double(temp);
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        return -charges[i] * charges[j] * Internal::quadraticInterpolate_diff1<ScalarType>(repDoubleSquareStep, erfcStep, x1, y[0], y[1], y[2], r);
    }

    template<class ScalarType, bool IsSmallCell>
    void RSpaceEwald<ScalarType, IsSmallCell>::makeTables() {
        for (size_t i = 2; i < erfc_table.getLength(); ++i) {
            const auto x = PlainScalar((i - 1) * ErfcTableStep);
            erfc_table[i] = erfc(x) / x * integralLimit;
        }
        erfc_table[0] = erfc_table[1] = erfc_table[2]; // Smooth out divergent erfc(0) / 0 
        erfcStep = PlainScalar(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        repDoubleSquareStep = reciprocal(square(erfcStep) * PlainScalar(2));
    }
    /**
     * Slow version functors are provided for debug use
     */
    template<class ScalarType, bool IsSmallCell>
    ScalarType RSpaceEwald<ScalarType, IsSmallCell>::pot_functor_slow(size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType x = r * integralLimit;
        return charges[i] * charges[j] * erfc(x) / r;
    }

    template<class ScalarType, bool IsSmallCell>
    ScalarType RSpaceEwald<ScalarType, IsSmallCell>::force_functor_slow(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        const ScalarType x = r * integralLimit;
        return charges[i] * charges[j] * (erfc(x) + x * exp(-square(x)) * PlainScalar(M_2_SQRTPI)) / r2;
    }

    template<class ScalarType, bool IsSmallCell>
    ScalarType RSpaceEwald<ScalarType, IsSmallCell>::rSpaceForceConstImpl1(ScalarType r) const {
        const ScalarType x = integralLimit * r;
        const ScalarType x2 = square(x);
        const ScalarType term1 = ScalarType(3) * erfc(x);
        const ScalarType term2 = ScalarType(2 * M_2_SQRTPI) * x * (ScalarType(1.5) + x2) / exp(x2);
        return (term1 + term2) / (square(square(r)) * r);
    }

    template<class ScalarType, bool IsSmallCell>
    ScalarType RSpaceEwald<ScalarType, IsSmallCell>::rSpaceForceConstImpl2(ScalarType r) const {
        const ScalarType x = integralLimit * r;
        const ScalarType term1 = erfc(x);
        const ScalarType term2 = ScalarType(M_2_SQRTPI) * x / exp(square(x));
        return -(term1 + term2) / (square(r) * r);
    }

    template<class ScalarType, bool IsSmallCell>
    typename RSpaceEwald<ScalarType, IsSmallCell>::BornChargeArray
    RSpaceEwald<ScalarType, IsSmallCell>::makeBornCharge(const Vector<ScalarType>& charges) {
        BornChargeArray result(charges.getLength(), Matrix3D(3, 3, ScalarType(0)));
        for (size_t i = 0; i < result.getLength(); ++i) {
            auto diag = result[i].diag();
            diag = charges[i];
        }
        return result;
    }
}

namespace Physica {
    template<class T, bool B>
    class Traits<Core::RSpaceEwald<T, B>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPotDependOnAtomIndex = true;
        constexpr static bool IsSmallCell = B;
    };
}
