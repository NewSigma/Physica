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
    template<class T, class ReprType, class VectorType>
    class MatrixVectorProduct<HubbardMatrix<T, ReprType>, VectorType>
            : public RValueVector<MatrixVectorProduct<HubbardMatrix<T, ReprType>, VectorType>> {
        using MatrixType = HubbardMatrix<T, ReprType>;
        using This = MatrixVectorProduct<MatrixType, VectorType>;
        using Base = RValueVector<This>;
        constexpr static unsigned int Dim = MatrixType::Dim;
        constexpr static unsigned int NumSite = MatrixType::NumSite;
        constexpr static unsigned int SiteDOF = MatrixType::SiteDOF;
        constexpr static bool IsTransInvariant = MatrixType::IsTransInvariant;
    public:
        using ScalarType = typename Base::ScalarType;
        using RealType = typename ScalarType::RealType;
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const RValueMatrix<MatrixType>& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(OtherDerived& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class T, class ReprType, class VectorType>
    template<class OtherDerived, class Executor>
    inline void MatrixVectorProduct<HubbardMatrix<T, ReprType>, VectorType>::assignTo(OtherDerived& target) const {
        static_assert(std::is_base_of<RValueVector<VectorType>, VectorType>::value, "[Error]: Invalid source type");
        static_assert(std::is_base_of<LValueVector<OtherDerived>, OtherDerived>::value, "[Error]: Invalid target type");
        assert(target.getLength() == getLength() && "[Error]: Dimensions do not match");

        target = RealType(0);
        if constexpr (std::is_same<Executor, ThreadExecutor>::value) {
            std::mutex mutex{};
            auto future = Executor::parallel_for([&](unsigned int thread) {
                const size_t length = getLength();
                Vector<ScalarType> local(length, 0);
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
                mat.dotImpl(target, vec.calc(i), i);
        }
    }

    template<class T, class ReprType, class VectorType>
    typename MatrixVectorProduct<HubbardMatrix<T, ReprType>, VectorType>::ScalarType
    MatrixVectorProduct<HubbardMatrix<T, ReprType>, VectorType>::calc(size_t index) const {
        static_assert(!IsTransInvariant && "[Error]: Not implemented");
        const ScalarType hop = -mat.getHoppingT();
        const auto state = mat.repr[index];
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
                    const size_t index1 = mat.repr[upOccupy1 ? state.hopUp(site, site1) : state.hopUp(site1, site)];
                    result += vec.calc(index1) * (upOccupy1 ? hopUp : -hopUp);
                }

                if (downOccupy1 != downOccupy2) {
                    const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                    const size_t index1 = mat.repr[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                    result += vec.calc(index1) * (downOccupy1 ? hopDown : -hopDown);
                }
            }, site);
            numRepel += upOccupy1 && downOccupy1;
        }
        result += vec.calc(index) * (mat.getRepelU() * RealType(numRepel));
        return result;
    }
}

namespace Physica {
    using namespace Core;

    template<class T, class ReprType, class VectorType>
    class Traits<MatrixVectorProduct<HubbardMatrix<T, ReprType>, VectorType>> {
        using MatrixType = HubbardMatrix<T, ReprType>;
        using T1 = typename VectorType::ScalarType;
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<T, T1>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
