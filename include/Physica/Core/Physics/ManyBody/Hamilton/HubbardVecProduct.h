/*
 * Copyright 2024 Weibo He.
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
    template<Scalar T, class ReprType, Vector U>
    class MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>
            : public RValueVector<MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>> {
        using MatrixType = HubbardMatrix<T, ReprType>;
        using This = MatrixVectorProduct<MatrixType, U>;
        using Base = RValueVector<This>;
        using FFTType = typename MatrixType::FFTType;
        using StateType = typename ReprType::StateType;
        constexpr static unsigned int Dim = MatrixType::Dim;
        constexpr static unsigned int NumSite = MatrixType::NumSite;
        constexpr static unsigned int SiteDOF = MatrixType::SiteDOF;
        constexpr static bool IsTransInvariant = MatrixType::IsTransInvariant;
    public:
        using ScalarType = typename Base::ScalarType;
        using RealType = typename ScalarType::RealType;
    private:
        const MatrixType& mat;
        const U& vec;
    public:
        MatrixVectorProduct(const MatrixType& mat_, const U& vec_) : mat(mat_), vec(vec_) {
            assert(mat.getCol() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const U& getRHS() const noexcept { return vec; }
    private:
        template<Vector V> void sumHopping(V& target, FFTType& fft, ScalarType factor, StateType psi) const;
        template<Vector V> void dotImpl(V& target, ScalarType factor, size_t index) const;
        /* Getters */
        [[nodiscard]] ScalarType getHoppingT() const noexcept { return mat.getHoppingT(); }
        [[nodiscard]] ScalarType getRepelU() const noexcept { return mat.getRepelU(); }
        [[nodiscard]] const ReprType& getRepr() const noexcept { return mat.getRepr(); }
    };

    template<Scalar T, class ReprType, Vector U>
    template<LVector V, class Executor>
    inline void MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>::assignTo(V& target) const {
        assert(target.getLength() == getLength() && "[Error]: Dimensions do not match");
        target = RealType(0);
        if constexpr (std::is_same<Executor, ThreadExecutor>::value) {
            std::mutex mutex{};
            auto future = Executor::parallel_for([&](unsigned int thread) {
                const size_t length = getLength();
                VectorND<ScalarType> local(length, 0);
                SparseVector<ScalarType> buffer(length, std::min(size_t(NumSite * SiteDOF), length));
                const auto range = Executor::splitJob(length, Executor::getNumThread(), thread);
                for (unsigned int i = range.first; i < range.second; ++i) {
                    mat.dotImpl(buffer, vec.calc(i), i);
                    local += buffer;
                    buffer.clear();
                }
                std::unique_lock locker(mutex);
                target += local;
            }, Executor::getNumThread());
            Executor::auto_wait(future);
        }
        else {
            const size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                dotImpl<V>(target, vec.calc(i), i);
        }
    }

    template<Scalar T, class ReprType, Vector U>
    typename MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>::ScalarType
    MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>::calc(size_t index) const {
        static_assert(!IsTransInvariant && "[Error]: Not implemented");
        const ScalarType hop = -mat.getHoppingT();
        const auto state = getRepr()[index];
        ScalarType result = 0;
        int numRepel = 0;
        for (int site = 0; site < int(NumSite); ++site) {
            const bool upOccupy1 = state.isUpOccupy(site);
            const bool downOccupy1 = state.isDownOccupy(site);
            mat.forNeighSites([this, &result, &state, hop, upOccupy1, downOccupy1](int site, int site1) {
                const bool upOccupy2 = state.isUpOccupy(site1);
                const bool downOccupy2 = state.isDownOccupy(site1);
                if (upOccupy1 != upOccupy2) {
                    const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                    const size_t index1 = getRepr()[upOccupy1 ? state.hopUp(site, site1) : state.hopUp(site1, site)];
                    result += vec.calc(index1) * (upOccupy1 ? hopUp : -hopUp);
                }

                if (downOccupy1 != downOccupy2) {
                    const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                    const size_t index1 = getRepr()[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                    result += vec.calc(index1) * (downOccupy1 ? hopDown : -hopDown);
                }
            }, site);
            numRepel += upOccupy1 && downOccupy1;
        }
        result += vec.calc(index) * (mat.getRepelU() * RealType(numRepel));
        return result;
    }

    template<Scalar T, class ReprType, Vector U>
    template<Vector V>
    void MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>::sumHopping(V& target, FFTType& fft, ScalarType factor, StateType psi) const {
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
        FFTType::transform(mat.planProvider, fft);
        const auto& repr = getRepr();
        const size_t index = repr[reducedPsi];
        target[index] += fft.getKSpace()[repr.getReducedK()] * sqrt(RealType(repr.getPeriods()[index])) * factor;
    }

    template<Scalar T, class ReprType, Vector U>
    template<Vector V>
    void MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>::dotImpl(V& target, ScalarType factor, size_t index) const {
        const auto state = getRepr()[index];
        int numRepel = 0;
        if constexpr (IsTransInvariant) {
            static_assert(Dim == 1 && "[Error]: Not implemented");
            const RealType normalizer = sqrt(RealType(getRepr().getPeriods()[index])) / RealType(NumSite);
            const ScalarType hop = -factor * normalizer * getHoppingT();

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
            const ScalarType hop = -factor * getHoppingT();
            for (int site = 0; site < int(NumSite); ++site) {
                const bool upOccupy1 = state.isUpOccupy(site);
                const bool downOccupy1 = state.isDownOccupy(site);
                mat.forNeighSites([this, &target, &state, hop, upOccupy1, downOccupy1](int site, int site1) {
                    const auto& repr = getRepr();
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
                }, site);
                numRepel += upOccupy1 && downOccupy1;
            }
        }
        target[index] += factor * (getRepelU() * RealType(numRepel));
    }
}

namespace Physica {
    template<Scalar T, class ReprType, Vector U>
    class Traits<MatrixVectorProduct<HubbardMatrix<T, ReprType>, U>> {
        using MatrixType = HubbardMatrix<T, ReprType>;
        using T1 = typename U::ScalarType;
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpRtnTy<T, T1>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
