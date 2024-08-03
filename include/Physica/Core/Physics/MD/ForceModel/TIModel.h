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

#include "Physica/Core/Physics/PhyConst.h"
#include "Hamonic.h"

namespace Physica::Core {
    /**
     * \class TIModel provides support to thermodynamic integration(TI).
     * 
     * Reference:
     * [1] Comput. Mater. Sci. 112, 333-341 (2016); https://doi.org/10.1016/j.commatsci.2015.10.050
     */
    template<class ForceModel>
    class TIModel {
        using ScalarType = typename Traits<ForceModel>::ScalarType;
        using HamonicType = Hamonic<ScalarType, 3>;
        using MDCellType = typename HamonicType::MDCellType;
        using PositionMatrix = typename HamonicType::PositionMatrix;

        ForceModel original;
        HamonicType hamonic;
        Vector<ScalarType> mass;
        ScalarType temperatureT;
        ScalarType refPotentialV;
        ScalarType refHelmholtzF;

        ScalarType lambda;
    public:
        TIModel(ForceModel original_, const MDCellType& refCell, Vector<ScalarType> springCoeffs, ScalarType temperatureT_);
        TIModel(const TIModel&) = default;
        TIModel(TIModel&&) noexcept = default;
        ~TIModel() = default;
        /* Operators */
        TIModel& operator=(TIModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType potentialV(const MDCellType& cell) const;
        [[nodiscard]] ScalarType deltaPotentialV(const MDCellType& cell) const;

        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result);

        void swap(TIModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getRefHelmholtzF() const noexcept { return refHelmholtzF; }
        /* Setters */
        void setLambda(ScalarType lambda_);
        /* Static members */
        [[nodiscard]] static Vector<ScalarType> makeSpringCoeffs(const Vector<ScalarType>& msd, ScalarType temperatureT);
    private:
        void updateRef();
    };

    template<class ForceModel>
    TIModel<ForceModel>::TIModel(ForceModel original_, const MDCellType& refCell, Vector<ScalarType> springCoeffs, ScalarType temperatureT_)
            : original(std::move(original_))
            , hamonic(refCell.getPos(), std::move(springCoeffs))
            , mass(refCell.getMassVec())
            , temperatureT(temperatureT_)
            , lambda(0) {
        updateRef();
        refPotentialV = original.potentialV(refCell);
    }
    
    template<class ForceModel>
    typename TIModel<ForceModel>::ScalarType TIModel<ForceModel>::potentialV(const MDCellType& cell) const {
        const ScalarType hamonicV = hamonic.potentialV(cell);
        if (lambda.isZero())
            return hamonicV;
        return lambda * original.potentialV(cell) + (ScalarType(1) - lambda) * hamonicV;
    }

    template<class ForceModel>
    typename TIModel<ForceModel>::ScalarType TIModel<ForceModel>::deltaPotentialV(const MDCellType& cell) const {
        return (original.potentialV(cell) - refPotentialV) - hamonic.potentialV(cell);
    }

    template<class ForceModel>
    template<class Executor>
    Vector<typename TIModel<ForceModel>::ScalarType> TIModel<ForceModel>::force(const MDCellType& cell) const {
        const Vector<ScalarType> hamonicF = hamonic.force(cell);
        if (lambda.isZero())
            return hamonicF;
        return lambda * original.template force<Executor>(cell) + (ScalarType(1) - lambda) * hamonicF;
    }

    template<class ForceModel>
    template<class VectorType, class Executor>
    void TIModel<ForceModel>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) {
        if (!lambda.isZero())
            original.template forceAsync<VectorType, Executor>(cell, result);
        const Vector<ScalarType> hamonicF = hamonic.template force<Executor>(cell);
        Executor::wait();
        result = lambda * result + (ScalarType(1) - lambda) * hamonicF;
    }

    template<class ForceModel>
    void TIModel<ForceModel>::swap(TIModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        original.swap(obj.original);
        hamonic.swap(obj.hamonic);
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
    Vector<typename TIModel<ForceModel>::ScalarType> TIModel<ForceModel>::makeSpringCoeffs(
            const Vector<ScalarType>& msd, ScalarType temperatureT) {
        const ScalarType factor = ScalarType(3 * PhyConst<AU>::boltzmannK) * temperatureT;
        return reciprocal(msd) * factor;
    }

    template<class ForceModel>
    void TIModel<ForceModel>::updateRef() {
        const ScalarType factor = ScalarType(PhyConst<AU>::reducedPlanck / PhyConst<AU>::boltzmannK) / temperatureT;
        const Vector<ScalarType> omegas = sqrt(divide(hamonic.getSpringCoeffs(), mass));
        refHelmholtzF = ScalarType(3 * PhyConst<AU>::boltzmannK) * temperatureT * ln(omegas * factor).sum();
    }
}

namespace Physica {
    template<class ForceModel>
    class Traits<Core::TIModel<ForceModel>> : public Traits<ForceModel> {};
}
