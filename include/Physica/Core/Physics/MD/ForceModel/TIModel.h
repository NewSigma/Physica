/*
 * Copyright 2024-2026 Weibo He.
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
#include "Harmonic.h"

namespace Physica {
    /**
     * \class TIModel provides support to thermodynamic integration(TI).
     * 
     * Reference:
     * [1] Comput. Mater. Sci. 112, 333-341 (2016); https://doi.org/10.1016/j.commatsci.2015.10.050
     */
    template<class ForceModel>
    class TIModel {
        using This = TIModel<ForceModel>;
        using ScalarType = Traits<ForceModel>::ScalarType;
        using HarmonicType = Harmonic<ScalarType, 3>;
        using MDCellType = HarmonicType::MDCellType;
        using PositionMatrix = HarmonicType::PositionMatrix;

        ForceModel original;
        HarmonicType harmonic;
        VectorND<ScalarType> mass;
        ScalarType temperatureT;
        ScalarType refPotentialV;
        ScalarType refHelmholtzF;

        ScalarType lambda;
    public:
        TIModel(ForceModel original_, const MDCellType& refCell, VectorND<ScalarType> springCoeffs, ScalarType temperatureT_);
        TIModel(const This&) = default;
        TIModel(This&&) noexcept = default;
        ~TIModel() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType potentialV(const MDCellType& cell) const;
        [[nodiscard]] ScalarType deltaPotentialV(const MDCellType& cell) const;

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<ScalarType> force(const MDCellType& cell) const;
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getRefHelmholtzF() const noexcept { return refHelmholtzF; }
        /* Setters */
        void setLambda(ScalarType lambda_);
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return ForceModel::isPeriodBoundary(); }
        [[nodiscard]] static VectorND<ScalarType> makeSpringCoeffs(const VectorND<ScalarType>& msd, ScalarType temperatureT);
    private:
        void updateRef();
    };

    template<class ForceModel>
    TIModel<ForceModel>::TIModel(ForceModel original_, const MDCellType& refCell, VectorND<ScalarType> springCoeffs, ScalarType temperatureT_)
            : original(std::move(original_))
            , harmonic(refCell.getPos(), std::move(springCoeffs))
            , mass(refCell.getMassVec())
            , temperatureT(temperatureT_)
            , refPotentialV(original.potentialV(refCell))
            , lambda(0) {
        updateRef();
    }
    
    template<class ForceModel>
    TIModel<ForceModel>::ScalarType TIModel<ForceModel>::potentialV(const MDCellType& cell) const {
        const ScalarType harmonicV = harmonic.potentialV(cell);
        if (lambda.isZero())
            return harmonicV;
        return lambda * original.potentialV(cell) + (ScalarType(1) - lambda) * harmonicV;
    }

    template<class ForceModel>
    TIModel<ForceModel>::ScalarType TIModel<ForceModel>::deltaPotentialV(const MDCellType& cell) const {
        return (original.potentialV(cell) - refPotentialV) - harmonic.potentialV(cell);
    }

    template<class ForceModel>
    template<ExecutePolicy P>
    VectorND<typename TIModel<ForceModel>::ScalarType> TIModel<ForceModel>::force(const MDCellType& cell) const {
        const VectorND<ScalarType> harmonicF = harmonic.force(cell);
        if (lambda.isZero())
            return harmonicF;
        return lambda * original.template force<P>(cell) + (ScalarType(1) - lambda) * harmonicF;
    }

    template<class ForceModel>
    template<ExecutePolicy P>
    void TIModel<ForceModel>::forceAsync(const MDCellType& cell, Vector auto& result) {
        if (!lambda.isZero())
            original.template forceAsync<P>(cell, result);
        const VectorND<ScalarType> harmonicF = harmonic.template force<P>(cell);
        if constexpr (P == GPU)
            Task<GPU>::wait();
        result = lambda * result.getDerived() + (ScalarType(1) - lambda) * harmonicF;
    }

    template<class ForceModel>
    void TIModel<ForceModel>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        original.swap(obj.original);
        harmonic.swap(obj.harmonic);
        mass.swap(obj.mass);
        temperatureT.swap(obj.temperatureT);
        refPotentialV.swap(obj.refPotentialV);
        refHelmholtzF.swap(obj.refHelmholtzF);
        lambda.swap(obj.lambda);
    }

    template<class ForceModel>
    void TIModel<ForceModel>::setLambda(ScalarType lambda_) {
        assert((ScalarType(0) <= lambda_) && (lambda_ <= ScalarType(1)) && "[Error]: Invalid value for thermodynamic integration");
        lambda = lambda_;
    }
    /**
     * Fast method to generate spring coefficients from MSD as [1] does.
     * 
     * Reference:
     * [1] Comput. Mater. Sci. 112, 333-341 (2016); https://doi.org/10.1016/j.commatsci.2015.10.050
     */
    template<class ForceModel>
    VectorND<typename TIModel<ForceModel>::ScalarType> TIModel<ForceModel>::makeSpringCoeffs(
            const VectorND<ScalarType>& msd, ScalarType temperatureT) {
        const ScalarType factor = ScalarType(3 * PhyConst<AU>::boltzmannK) * temperatureT;
        return reciprocal(msd) * factor;
    }

    template<class ForceModel>
    void TIModel<ForceModel>::updateRef() {
        const ScalarType factor = ScalarType(PhyConst<AU>::reducedPlanck / PhyConst<AU>::boltzmannK) / temperatureT;
        const VectorND<ScalarType> omegas = sqrt(divide(harmonic.getSpringCoeffs(), mass));
        refHelmholtzF = ScalarType(3 * PhyConst<AU>::boltzmannK) * temperatureT * ln(omegas * factor).sum();
    }
}

namespace Physica {
    template<class ForceModel>
    class Traits<TIModel<ForceModel>> : public Traits<ForceModel> {};
}
