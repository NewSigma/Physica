/*
 * Copyright 2023 WeiBo He.
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

#include "RSpaceEwald.h"
#include "Physica/Core/Physics/MD/ForceModel/PairModel.cuh"

namespace Physica::Core {
    template<class ScalarType, bool IsSmallCell>
    class device_obj<RSpaceEwald<ScalarType, IsSmallCell>> : public device_obj<PairModel<RSpaceEwald<ScalarType, IsSmallCell>>> {
        static_assert(!IsSmallCell, "[Error]: Small cell does not apply to ewald because self interaction");
        using host_obj = RSpaceEwald<ScalarType, IsSmallCell>;
        using This = device_obj<host_obj>;
        using Base = device_obj<PairModel<host_obj>>;
        using DeviceVector = device_obj<Vector<ScalarType>>;
    public:
        using Base::Dim;
        using typename Base::PlainScalar;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::PositionMatrix;
        using DeviceLattice = device_obj<LatticeMatrix>;
        using SearchRangeType = typename device_obj<PeriodicCell<ScalarType, Dim>>::SearchRangeType;
        constexpr static size_t ErfcTableSize = host_obj::ErfcTableSize;
        constexpr static double ErfcTableStep = host_obj::ErfcTableStep;
        constexpr static double SumPrec = host_obj::SumPrec;
    private:
        DeviceLattice lattice;
        DeviceLattice repLatt;
        InvLatticeMatrix invLatt;
        DeviceVector charges;
        DeviceVector erfc_table;
        ScalarType volume;
        ScalarType inv_volume;
        ScalarType integralLimit;
        ScalarType erfcStep;
        ScalarType repErfcStep;
        ScalarType repDoubleSquareStep;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
    public:
        device_obj() = default;
        device_obj(const LatticeMatrix& lattice_, const Vector<ScalarType>& charges_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operator */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Executor>
        [[nodiscard]] inline Vector<ScalarType> force_short(const PositionMatrix& pos);

        [[nodiscard]] inline LatticeMatrix virial(const PositionMatrix& pos);
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return charges.getLength(); }
        [[nodiscard]] __host__ __device__ ScalarType getVolume() const noexcept { return volume; }
        [[nodiscard]] __host__ __device__ ScalarType getInvVolume() const noexcept { return inv_volume; }
        [[nodiscard]] __host__ __device__ ScalarType getIntegralLimit() const noexcept { return integralLimit; }
        [[nodiscard]] __host__ __device__ ScalarType getRSpaceCutoff() const noexcept { return Base::getCutoff(); }
        [[nodiscard]] __host__ __device__ ScalarType getSquaredRSpaceCutoff() const noexcept { return Base::getSquaredCutoff(); }
        /* Setters */
        void setLattice(const LatticeMatrix& lattice_);
        void setIntegralLimit(ScalarType integralLimit_);
    protected:
        [[nodiscard]] inline ScalarType calcSelfE() const;
        [[nodiscard]] inline ScalarType calcGammaPointE() const;
        [[nodiscard]] __device__ inline ScalarType pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
        [[nodiscard]] __device__ inline ScalarType force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const;
    private:
        /* Operations */
        void makeTables();
        using Base::potentialEnergy;
        using Base::force_short;
        using Base::virial;
        /* Getters */
        using Base::getCutoff;
        using Base::getSquaredCutoff;
        /* Friends */
        friend class device_obj<PairModel<host_obj>>;
    };

    template<class ScalarType, bool IsSmallCell>
    device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::device_obj(const LatticeMatrix& lattice_, const Vector<ScalarType>& charges_)
            : Base(charges_.getLength()), charges(charges_), erfc_table(ErfcTableSize + 1) {
        setLattice(lattice_);
    }

    template<class ScalarType, bool IsSmallCell>
    template<class Executor>
    inline Vector<ScalarType> device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::force_short(const PositionMatrix& pos) {
        static_assert(std::is_same<Executor, CudaExecutor>::value, "[Error]: Invalid executor");
        static_assert(!IsSmallCell, "[Error]: Small cell does not apply to ewald because self interaction");
        const Vector<ScalarType> rSpaceSum = Base::template force<CudaExecutor>(lattice.toHost(), invLatt, pos);
        return rSpaceSum;
    }

    template<class ScalarType, bool IsSmallCell>
    inline typename device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::LatticeMatrix
    device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::virial(const PositionMatrix& pos) {
        return Base::virial(lattice.toHost(), invLatt, pos);
    }

    template<class ScalarType, bool IsSmallCell>
    void device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        lattice.swap(obj.lattice);
        repLatt.swap(obj.repLatt);
        invLatt.swap(obj.invLatt);
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
    void device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::setLattice(const LatticeMatrix& lattice_) {
        assert(charges.getLength() != 0 && "[Error]: Charges should be initialized before lattice update");
        lattice = lattice_;
        repLatt = PeriodicCell<ScalarType, Dim>::makeRepLattice(lattice_);
        invLatt = lattice_.inverse();
        volume = PeriodicCell<ScalarType, Dim>::getVolume(lattice_);
        inv_volume = reciprocal(volume);

        const ScalarType averageCellSize = cbrt(ScalarType(volume));
        const ScalarType estimate = sqrt(cbrt(ScalarType(getNumParticle())) * ScalarType(M_PI)) / averageCellSize;
        setIntegralLimit(estimate);
    }

    template<class ScalarType, bool IsSmallCell>
    void device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::setIntegralLimit(ScalarType integralLimit_) {
        const auto hostRepLatt = repLatt.toHost();
        const ScalarType heightX_2Pi = reciprocal(hostRepLatt.row(0).norm());
        const ScalarType heightY_2Pi = reciprocal(hostRepLatt.row(1).norm());
        const ScalarType heightZ_2Pi = reciprocal(hostRepLatt.row(2).norm());
        constexpr double factor1 = 2 * M_PI * (1 - std::numeric_limits<ScalarType>::epsilon()); //To avoid rSpaceCutoff larger than max value
        const ScalarType maxRSpaceCutoff = std::min(heightX_2Pi, std::min(heightY_2Pi, heightZ_2Pi)) * PlainScalar(factor1);
        const ScalarType minLimit = ScalarType(SumPrec) / maxRSpaceCutoff;
        integralLimit = std::max(integralLimit_, minLimit).getValue();

        const PlainScalar rSpaceCutoff = PlainScalar(SumPrec) / integralLimit.getValue();
        rSpaceSumRange = PeriodicCell<ScalarType, Dim>::estimateRange(lattice.toHost(), rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<ScalarType, Dim>::estimateRange(hostRepLatt, PlainScalar(SumPrec * 2) * integralLimit.getValue());
        makeTables();
        Base::setCutoff(rSpaceCutoff);
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::calcSelfE() const {
        return square(charges.toHost()).sum() * integralLimit / sqrt(PlainScalar(M_PI));
    }

    template<class ScalarType, bool IsSmallCell>
    inline ScalarType device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::calcGammaPointE() const {
        return square(charges.toHost().sum()) * PlainScalar(-M_PI) / (PlainScalar(2) * square(integralLimit)) * inv_volume;
    }
    /**
     * Optimize: make use of x1, x2, x3 are equal distance
     */
    template<class ScalarType, bool IsSmallCell>
    __device__ inline ScalarType device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::pot_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = temp.getTrivial();
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        const ScalarType interp = Internal::quadraticInterpolate<ScalarType>(x1 - erfcStep, x1, x1 + erfcStep, y[0], y[1], y[2], r);
        return charges[i] * charges[j] * interp;
    }

    template<class ScalarType, bool IsSmallCell>
    __device__ inline ScalarType device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::force_functor(
            size_t i, size_t j, ScalarType r, [[maybe_unused]] ScalarType r2) const {
        const ScalarType temp = r * repErfcStep + PlainScalar(0.5);
        const int index = temp.getTrivial();
        const ScalarType x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        return -charges[i] * charges[j] * Internal::quadraticInterpolate_diff1<ScalarType>(repDoubleSquareStep, erfcStep, x1, y[0], y[1], y[2], r);
    }

    template<class ScalarType, bool IsSmallCell>
    void device_obj<RSpaceEwald<ScalarType, IsSmallCell>>::makeTables() {
        Vector<ScalarType> hostErfcTable(erfc_table.getLength());
        for (size_t i = 2; i < hostErfcTable.getLength(); ++i) {
            const auto x = PlainScalar((i - 1) * ErfcTableStep);
            hostErfcTable[i] = erfc(x) / x * integralLimit;
        }
        hostErfcTable[0] = hostErfcTable[1] = hostErfcTable[2]; // Smooth out divergent erfc(0) / 0
        hostErfcTable.toDevice(erfc_table);
        erfcStep = PlainScalar(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        repDoubleSquareStep = reciprocal(square(erfcStep) * PlainScalar(2));
    }
}
