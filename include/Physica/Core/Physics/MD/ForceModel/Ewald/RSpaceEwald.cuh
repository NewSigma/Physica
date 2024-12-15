/*
 * Copyright 2023-2024 Weibo He.
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
    template<Scalar T, bool IsSmallCell>
    class device_obj<RSpaceEwald<T, IsSmallCell>> : public device_obj<PairModel<RSpaceEwald<T, IsSmallCell>>> {
        static_assert(!IsSmallCell, "[Error]: Small cell does not apply to ewald because self interaction");
        using host_obj = RSpaceEwald<T, IsSmallCell>;
        using This = device_obj<host_obj>;
        using Base = device_obj<PairModel<host_obj>>;
        using DeviceVector = device_obj<VectorND<T>>;
    public:
        using Base::Dim;
        using typename Base::ValueType;
        using typename Base::LatticeMatrix;
        using typename Base::InvLatticeMatrix;
        using typename Base::PositionMatrix;
        using DeviceLattice = device_obj<LatticeMatrix>;
        using SearchRangeType = device_obj<PeriodicCell<T, Dim>>::SearchRangeType;
        using BornChargeArray = host_obj::BornChargeArray;
        constexpr static size_t ErfcTableSize = host_obj::ErfcTableSize;
        constexpr static double ErfcTableStep = host_obj::ErfcTableStep;
        constexpr static double SumPrec = host_obj::SumPrec;
    private:
        DeviceLattice lattice;
        DeviceLattice repLatt;
        InvLatticeMatrix invLatt;
        DeviceVector charges;
        DeviceVector erfc_table;
        T volume;
        T inv_volume;
        T integralLimit;
        T erfcStep;
        T repErfcStep;
        T repDoubleSquareStep;
        SearchRangeType rSpaceSumRange;
        SearchRangeType kSpaceSumRange;
    public:
        device_obj() = default;
        device_obj(const LatticeMatrix& lattice_, const VectorND<T>& charges_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operator */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Executor>
        [[nodiscard]] inline VectorND<T> force_short(const PositionMatrix& pos);

        [[nodiscard]] inline LatticeMatrix virial(const PositionMatrix& pos);
        [[nodiscard]] BornChargeArray calcBornCharge() const { return host_obj::makeBornCharge(charges); }
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const DeviceLattice& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const DeviceLattice& getRepLattice() const noexcept { return repLatt; }
        [[nodiscard]] const DeviceVector& getCharges() const noexcept { return charges; }
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return charges.getLength(); }
        [[nodiscard]] __host__ __device__ T getVolume() const noexcept { return volume; }
        [[nodiscard]] __host__ __device__ T getInvVolume() const noexcept { return inv_volume; }
        [[nodiscard]] __host__ __device__ T getIntegralLimit() const noexcept { return integralLimit; }
        [[nodiscard]] __host__ __device__ T getRSpaceCutoff() const noexcept { return Base::getCutoff(); }
        [[nodiscard]] __host__ __device__ T getSquaredRSpaceCutoff() const noexcept { return Base::getSquaredCutoff(); }
        [[nodiscard]] const SearchRangeType& getRSpaceSumRange() const noexcept { return rSpaceSumRange; }
        [[nodiscard]] const SearchRangeType& getKSpaceSumRange() const noexcept { return kSpaceSumRange; }
        /* Setters */
        void setLattice(const LatticeMatrix& lattice_);
        void setIntegralLimit(T integralLimit_);
    protected:
        [[nodiscard]] inline T calcSelfE() const;
        [[nodiscard]] inline T calcGammaPointE() const;
        [[nodiscard]] __device__ inline T pot_functor(size_t i, size_t j, T r, T r2) const;
        [[nodiscard]] __device__ inline T force_functor(size_t i, size_t j, T r, T r2) const;
    private:
        /* Operations */
        void makeTables();
        using Base::potentialV;
        using Base::virial;
        /* Getters */
        using Base::getCutoff;
        using Base::getSquaredCutoff;
        /* Friends */
        friend class device_obj<PairModel<host_obj>>;
    };

    template<Scalar T, bool IsSmallCell>
    device_obj<RSpaceEwald<T, IsSmallCell>>::device_obj(const LatticeMatrix& lattice_, const VectorND<T>& charges_)
            : Base(charges_.getLength()), charges(charges_), erfc_table(ErfcTableSize + 1) {
        setLattice(lattice_);
    }

    template<Scalar T, bool IsSmallCell>
    template<class Executor>
    inline VectorND<T> device_obj<RSpaceEwald<T, IsSmallCell>>::force_short(const PositionMatrix& pos) {
        static_assert(std::is_same<Executor, CUDAExecutor>::value, "[Error]: Invalid executor");
        static_assert(!IsSmallCell, "[Error]: Small cell does not apply to ewald because self interaction");
        const VectorND<T> rSpaceSum = Base::template force<CUDAExecutor>(lattice.toHost(), invLatt, pos);
        return rSpaceSum;
    }

    template<Scalar T, bool IsSmallCell>
    inline device_obj<RSpaceEwald<T, IsSmallCell>>::LatticeMatrix
    device_obj<RSpaceEwald<T, IsSmallCell>>::virial(const PositionMatrix& pos) {
        return Base::virial(lattice.toHost(), invLatt, pos);
    }

    template<Scalar T, bool IsSmallCell>
    void device_obj<RSpaceEwald<T, IsSmallCell>>::swap(device_obj& __restrict obj) noexcept {
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

    template<Scalar T, bool IsSmallCell>
    void device_obj<RSpaceEwald<T, IsSmallCell>>::setLattice(const LatticeMatrix& lattice_) {
        assert(charges.getLength() != 0 && "[Error]: Charges should be initialized before lattice update");
        lattice = lattice_;
        repLatt = PeriodicCell<T, Dim>::makeRepLattice(lattice_);
        invLatt = lattice_.inverse();
        volume = PeriodicCell<T, Dim>::getVolume(lattice_);
        inv_volume = reciprocal(volume);

        const T averageCellSize = cbrt(T(volume));
        const T estimate = sqrt(cbrt(T(getNumParticle())) * T(M_PI)) / averageCellSize;
        setIntegralLimit(estimate);
    }

    template<Scalar T, bool IsSmallCell>
    void device_obj<RSpaceEwald<T, IsSmallCell>>::setIntegralLimit(T integralLimit_) {
        const auto hostRepLatt = repLatt.toHost();
        const T heightX_2Pi = reciprocal(hostRepLatt.row(0).norm());
        const T heightY_2Pi = reciprocal(hostRepLatt.row(1).norm());
        const T heightZ_2Pi = reciprocal(hostRepLatt.row(2).norm());
        constexpr double factor1 = 2 * M_PI * (1 - std::numeric_limits<T>::epsilon()); //To avoid rSpaceCutoff larger than max value
        const T maxRSpaceCutoff = std::min(heightX_2Pi, std::min(heightY_2Pi, heightZ_2Pi)) * ValueType(factor1);
        const T minLimit = T(SumPrec) / maxRSpaceCutoff;
        integralLimit = std::max(integralLimit_, minLimit).value();

        const ValueType rSpaceCutoff = ValueType(SumPrec) / integralLimit.value();
        rSpaceSumRange = PeriodicCell<T, Dim>::estimateRange(lattice.toHost(), rSpaceCutoff);
        kSpaceSumRange = PeriodicCell<T, Dim>::estimateRange(hostRepLatt, ValueType(SumPrec * 2) * integralLimit.value());
        makeTables();
        Base::setCutoff(rSpaceCutoff);
    }

    template<Scalar T, bool IsSmallCell>
    inline T device_obj<RSpaceEwald<T, IsSmallCell>>::calcSelfE() const {
        return square(charges.toHost()).sum() * integralLimit / sqrt(ValueType(M_PI));
    }

    template<Scalar T, bool IsSmallCell>
    inline T device_obj<RSpaceEwald<T, IsSmallCell>>::calcGammaPointE() const {
        return square(charges.toHost().sum()) * ValueType(-M_PI) / (ValueType(2) * square(integralLimit)) * inv_volume;
    }
    /**
     * Optimize: make use of x1, x2, x3 are equal distance
     */
    template<Scalar T, bool IsSmallCell>
    __device__ inline T device_obj<RSpaceEwald<T, IsSmallCell>>::pot_functor(
            size_t i, size_t j, T r, [[maybe_unused]] T r2) const {
        const T temp = r * repErfcStep + ValueType(0.5);
        const int index = temp.toMachine();
        const T x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        const T interp = Internal::quadraticInterpolate<T>(x1 - erfcStep, x1, x1 + erfcStep, y[0], y[1], y[2], r);
        return charges[i] * charges[j] * interp;
    }

    template<Scalar T, bool IsSmallCell>
    __device__ inline T device_obj<RSpaceEwald<T, IsSmallCell>>::force_functor(
            size_t i, size_t j, T r, [[maybe_unused]] T r2) const {
        const T temp = r * repErfcStep + ValueType(0.5);
        const int index = temp.toMachine();
        const T x1 = erfcStep * floor(temp);
        auto y = erfc_table.template segment<3>(index, index + 3);
        return -charges[i] * charges[j] * Internal::quadraticInterpolate_diff1<T>(repDoubleSquareStep, erfcStep, x1, y[0], y[1], y[2], r);
    }

    template<Scalar T, bool IsSmallCell>
    void device_obj<RSpaceEwald<T, IsSmallCell>>::makeTables() {
        VectorND<T> hostErfcTable(erfc_table.getLength());
        for (size_t i = 2; i < hostErfcTable.getLength(); ++i) {
            const auto x = ValueType((i - 1) * ErfcTableStep);
            hostErfcTable[i] = erfc(x) / x * integralLimit;
        }
        hostErfcTable[0] = hostErfcTable[1] = hostErfcTable[2]; // Smooth out divergent erfc(0) / 0
        hostErfcTable.toDevice(erfc_table);
        erfcStep = ValueType(ErfcTableStep) / integralLimit;
        repErfcStep = reciprocal(erfcStep);
        repDoubleSquareStep = reciprocal(square(erfcStep) * ValueType(2));
    }
}

namespace Physica {
    template<Scalar T, bool IsSmallCell>
    class Traits<Core::device_obj<RSpaceEwald<T, IsSmallCell>>> : public Traits<RSpaceEwald<T, IsSmallCell>> {};
}
