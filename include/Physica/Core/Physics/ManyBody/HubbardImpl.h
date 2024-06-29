/*
 * Copyright 2024 WeiBo He.
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
    template<class ScalarType, class ReprType>
    Hubbard<ScalarType, ReprType>::Hubbard(
            DimArray superSize, unsigned int numUnitCellSite, ReprType repr_, RealType hoppingT_, RealType repelU_)
            : Base(std::move(superSize), numUnitCellSite)
            , repr(std::move(repr_))
            , hoppingT(hoppingT_)
            , repelU(repelU_)
            , planProvider(NumSite, PlanFlag::Estimate) {
        if constexpr (Dim > 1)
            hoppingMatrix = makeHoppingMatrix();
    }

    template<class ScalarType, class ReprType>
    template<class SourceVector, class TargetVector, class Executor>
    void Hubbard<ScalarType, ReprType>::dot(const SourceVector& source, TargetVector& target) const {
        static_assert(std::is_base_of<RValueVector<SourceVector>, SourceVector>::value, "[Error]: Invalid source type");
        static_assert(std::is_base_of<LValueVector<TargetVector>, TargetVector>::value, "[Error]: Invalid target type");
        const size_t length = Base::getColumn();
        assert(source.getLength() == length && "[Error]: Dimensions do not match");
        assert(source.getLength() == target.getLength() && "[Error]: Dimensions do not match");

        target = RealType(0);
        if constexpr (std::is_same<Executor, ThreadExecutor>::value) {
            std::mutex mutex{};
            auto future = Executor::parallel_for([&, length](unsigned int thread) {
                VectorType local(length, 0);
                SparseType buffer(length, std::min(size_t(NumSite * SiteDOF), length));
                const auto range = Executor::splitJob(length, Executor::getNumThread(), thread);
                for (unsigned int i = range.first; i < range.second; ++i) {
                    if constexpr (Dim == 1)
                        dotImpl1D(buffer, source.calc(i), i);
                    else
                        dotImplND(buffer, source.calc(i), i);
                    local += buffer;
                    buffer.clear();
                }
                std::unique_lock locker(mutex);
                target += local;
            }, Executor::getNumThread());
            Executor::auto_wait(future);
        }
        else {
            for (size_t i = 0; i < length; ++i) {
                if constexpr (Dim == 1)
                    dotImpl1D(target, source.calc(i), i);
                else
                    dotImplND(target, source.calc(i), i);
            }
        }
    }

    template<class ScalarType, class ReprType>
    ScalarType Hubbard<ScalarType, ReprType>::calc(size_t row, size_t col) const {
        if constexpr (IsTransInvariant) {
            const auto& periods = repr.getPeriods();
            const auto psi1 = repr[row];
            if (row == col) {
                const bool flag = (repr.getKIndex() == 0) || (periods[row] == NumSite);
                return flag ? repelElem(psi1) : RealType(0);
            }

            auto fft = FFTType::makeEmptyFFT(NumSite);
            {
                auto& rSpace = fft.getRSpace();
                auto psi2 = repr[col];
                int sign = 1;
                for (int i = 0; i < int(NumSite); ++i) {
                    RealType elem = 0;
                    if (psi1 != psi2)
                        elem = hoppingElem(psi1, psi2) * RealType(sign);
                    rSpace[i] = elem;

                    sign *= psi2.lShiftSign();
                    psi2 <<= 1;
                }
            }
            FFTType::transform(planProvider, fft);
            const RealType normalizer = sqrt(RealType(periods[row] * periods[col])) / RealType(NumSite);
            return fft.getKSpace()[repr.getReducedK()] * normalizer;
        }
        else {
            const auto colState = repr[col];
            if (row == col)
                return repelElem(colState);
            return hoppingElem(repr[row], colState);
        }
    }

    template<class ScalarType, class ReprType>
    void Hubbard<ScalarType, ReprType>::swap(Hubbard& __restrict obj) noexcept {
        Base::swap(obj);
        repr.swap(obj.repr);
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        planProvider.swap(obj.planProvider);
    }

    template<class ScalarType, class ReprType>
    inline typename Hubbard<ScalarType, ReprType>::RealType Hubbard<ScalarType, ReprType>::repelElem(StateType psi) const {
        return repelU * RealType(psi.getNumPairedParticle());
    }

    template<class ScalarType, class ReprType>
    typename Hubbard<ScalarType, ReprType>::RealType
    Hubbard<ScalarType, ReprType>::hoppingElem(StateType rowPsi, StateType colPsi) const {
        int count = 0;
        if constexpr (Dim == 1) {
            for (int site = 0; site < int(NumSite); ++site) {
                const auto site1 = (site + 1) % NumSite;
                const int signUp = colPsi.hopUpSign(site, site1);
                const int signDown = colPsi.hopDownSign(site, site1);
                int n1 = 0;
                n1 += rowPsi == colPsi.hopUp(site, site1);
                n1 -= rowPsi == colPsi.hopUp(site1, site);
                int n2 = 0;
                n2 += rowPsi == colPsi.hopDown(site, site1);
                n2 -= rowPsi == colPsi.hopDown(site1, site);
                count += n1 * signUp + n2 * signDown;
            }
        }
        else {
            Base::forSiteInLattice([this, &count, rowPsi, colPsi](IndexType index) {
                const auto& dims = Base::getDims();
                const int site = IndexType::toIndex1D(dims, index);
                for (unsigned int dim = 0; dim < Dim; ++dim) {
                    IndexType index1 = index;
                    index1[dim] = (index1[dim] + 1) % Base::getSuperSize()[dim];

                    const int site1 = IndexType::toIndex1D(dims, index1);
                    const int signUp = colPsi.hopUpSign(site, site1);
                    const int signDown = colPsi.hopDownSign(site, site1);
                    int n1 = 0;
                    n1 += rowPsi == colPsi.hopUp(site, site1);
                    n1 -= rowPsi == colPsi.hopUp(site1, site);
                    int n2 = 0;
                    n2 += rowPsi == colPsi.hopDown(site, site1);
                    n2 -= rowPsi == colPsi.hopDown(site1, site);
                    count += n1 * signUp + n2 * signDown;
                }
            });
        }
        return RealType(-count) * hoppingT;
    }

    template<class ScalarType, class ReprType>
    void Hubbard<ScalarType, ReprType>::sumHopping(SparseType& buffer, FFTType& fft, ScalarType factor, StateType psi) const {
        if (psi.isVacuum())
            return;
        const auto reducedPsi = psi.transReduce();
        auto& rSpace = fft.getRSpace();
        int sign = 1;
        for (int i = 0; i < int(NumSite); ++i) {
            rSpace[i] = RealType(reducedPsi == psi ? sign : 0);
            sign *= psi.lShiftSign();
            psi <<= 1;
        }
        FFTType::transform(planProvider, fft);
        const size_t index = repr[reducedPsi];
        buffer[index] += fft.getKSpace()[repr.getReducedK()] * sqrt(RealType(repr.getPeriods()[index])) * factor;
    }

    template<class ScalarType, class ReprType>
    typename Hubbard<ScalarType, ReprType>::HoppingMatrix Hubbard<ScalarType, ReprType>::makeHoppingMatrix() {
        HoppingMatrix result(NumSite);
        Base::forSiteInLattice([this, &result](IndexType index) {
            const auto& dims = Base::getDims();
            Utils::Array<unsigned char> hopTargets{};
            hopTargets.reserve(NumSite * Dim * 2);
            for (unsigned int dim = 0; dim < Dim; ++dim) {
                IndexType index1 = index;
                index1[dim] = (index1[dim] + 1) % Base::getSuperSize()[dim];
                hopTargets.append(IndexType::toIndex1D(dims, index1));
            }
            hopTargets.squeeze();
            const auto site = IndexType::toIndex1D(dims, index);
            result[site] = std::move(hopTargets);
        });
        return result;
    }

    template<class ScalarType, class ReprType>
    template<class TargetVector>
    void Hubbard<ScalarType, ReprType>::dotImpl1D(TargetVector& target, ScalarType factor, size_t index) const {
        ScalarType hop = -factor * hoppingT;
        if constexpr (IsTransInvariant) {
            const RealType normalizer = sqrt(RealType(repr.getPeriods()[index])) / RealType(NumSite);
            hop *= normalizer;
        }

        const auto state = repr[index];
        int numRepel = 0;
        if constexpr (IsTransInvariant) {
            auto fft = FFTType::makeEmptyFFT(NumSite);
            for (int site = 0; site < int(NumSite); ++site) {
                const auto site1 = (site + 1) % NumSite;
                const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                sumHopping(target, fft, hopUp, state.hopUp(site, site1));
                sumHopping(target, fft, -hopUp, state.hopUp(site1, site));
                sumHopping(target, fft, hopDown, state.hopDown(site, site1));
                sumHopping(target, fft, -hopDown, state.hopDown(site1, site));
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            }
        }
        else {
            bool upOccupy1 = state.isUpOccupy(0);
            bool downOccupy1 = state.isDownOccupy(0);
            for (int site = 0; site < int(NumSite); ++site) {
                const auto site1 = (site + 1) % NumSite;
                const bool upOccupy2 = state.isUpOccupy(site1);
                const bool downOccupy2 = state.isDownOccupy(site1);
                if (upOccupy1 != upOccupy2) {
                    const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                    const size_t index = repr[upOccupy1 ? state.hopUp(site, site1) : state.hopUp(site1, site)];
                    target[index] += upOccupy1 ? hopUp : -hopUp;
                }

                if (downOccupy1 != downOccupy2) {
                    const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                    const size_t index = repr[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                    target[index] += downOccupy1 ? hopDown : -hopDown;
                }
                numRepel += upOccupy1 && downOccupy1;
                upOccupy1 = upOccupy2;
                downOccupy1 = downOccupy2;
            }
        }
        target[index] += factor * (repelU * RealType(numRepel));
    }

    template<class ScalarType, class ReprType>
    template<class TargetVector>
    void Hubbard<ScalarType, ReprType>::dotImplND(TargetVector& target, ScalarType factor, size_t index) const {
        static_assert(!IsTransInvariant && "[Error]: Not implemented");
        const ScalarType hop = -factor * hoppingT;
        const auto state = repr[index];
        int numRepel = 0;
        for (int site = 0; site < int(NumSite); ++site) {
            const auto& hopTargets = hoppingMatrix[site];
            const bool upOccupy1 = state.isUpOccupy(site);
            const bool downOccupy1 = state.isDownOccupy(site);
            for (int site1 : hopTargets) {
                const bool upOccupy2 = state.isUpOccupy(site1);
                const bool downOccupy2 = state.isDownOccupy(site1);
                if (upOccupy1 != upOccupy2) {
                    const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                    const size_t index = repr[upOccupy1 ? state.hopUp(site, site1) : state.hopUp(site1, site)];
                    target[index] += upOccupy1 ? hopUp : -hopUp;
                }

                if (downOccupy1 != downOccupy2) {
                    const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                    const size_t index = repr[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                    target[index] += downOccupy1 ? hopDown : -hopDown;
                }
            }
            numRepel += upOccupy1 && downOccupy1;
        }
        target[index] += factor * (repelU * RealType(numRepel));
    }
}
