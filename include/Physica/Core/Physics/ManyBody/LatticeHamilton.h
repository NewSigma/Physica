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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>

namespace Physica::Core {
    template<class Derived>
    class LatticeHamilton : public RValueMatrix<Derived> {
        using This = LatticeHamilton<Derived>;
        using Base = RValueMatrix<Derived>;
        using DerivedTraits = Traits<Derived>;
    public:
        using ScalarType = typename DerivedTraits::ScalarType;
        using ReprType = typename DerivedTraits::ReprType;
        using StateType = typename ReprType::StateType;
        using IndexType = typename StateType::IndexType;

        constexpr static unsigned int Dim = ReprType::Dim;
        constexpr static unsigned int NumSite = StateType::NumSite;
        using DimArray = Utils::Array<unsigned int, Dim>;
    private:
        DimArray superSize;
        unsigned char numUnitCellSite;
    public:
        ~LatticeHamilton() = default;
        /* Operations */
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        void swap(LatticeHamilton& __restrict obj) noexcept;

        template<class Functor> void forSiteInLattice(Functor func) const;
        /* Getters */
        [[nodiscard]] const DimArray& getSuperSize() const noexcept { return superSize; }
        [[nodiscard]] unsigned char getNumUnitCellSite() const noexcept { return numUnitCellSite; }
        [[nodiscard]] IndexType getDims() const noexcept;
        [[nodiscard]] const ReprType& getRepr() const noexcept { return Base::getDerived().getRepr(); }
        [[nodiscard]] size_t getNumState() const noexcept { return Base::getDerived().getNumState(); }
        [[nodiscard]] size_t getRow() const noexcept { return getNumState(); }
        [[nodiscard]] size_t getColumn() const noexcept { return getNumState(); }
    protected:
        LatticeHamilton() = default;
        LatticeHamilton(DimArray superSize_, unsigned char numSitePerCell_);
        LatticeHamilton(const LatticeHamilton&) = default;
        LatticeHamilton(LatticeHamilton&&) noexcept = default;
        /* Operators */
        LatticeHamilton& operator=(LatticeHamilton obj) noexcept { swap(obj); return *this; }
    private:
        [[nodiscard]] bool checkNumSite() const noexcept;
    };

    template<class Derived>
    LatticeHamilton<Derived>::LatticeHamilton(DimArray superSize_, unsigned char numSitePerCell_)
            : superSize(std::move(superSize_)), numUnitCellSite(numSitePerCell_) {
        assert(checkNumSite() && "[Error]: Inconsistent site number");
    }

    template<class Derived>
    void LatticeHamilton<Derived>::swap(LatticeHamilton& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        superSize.swap(obj.superSize);
        std::swap(numUnitCellSite, obj.numUnitCellSite);
    }

    template<class Derived>
    template<class Functor>
    void LatticeHamilton<Derived>::forSiteInLattice(Functor func) const {
        if constexpr (Dim == 1) {
            for (unsigned char x = 0; x < superSize[0]; ++x)
                for (unsigned char site = 0; site < numUnitCellSite; ++site)
                    func(IndexType{x, site});
        }
        if constexpr (Dim == 2) {
            for (unsigned char x = 0; x < superSize[0]; ++x)
                for (unsigned char y = 0; y < superSize[1]; ++y)
                    for (unsigned char site = 0; site < numUnitCellSite; ++site)
                        func(IndexType{x, y, site});
        }
        if constexpr (Dim == 3) {
            for (unsigned char x = 0; x < superSize[0]; ++x)
                for (unsigned char y = 0; y < superSize[1]; ++y)
                    for (unsigned char z = 0; z < superSize[2]; ++z)
                        for (unsigned char site = 0; site < numUnitCellSite; ++site)
                            func(IndexType{x, y, z, site});
        }
    }

    template<class Derived>
    typename LatticeHamilton<Derived>::IndexType LatticeHamilton<Derived>::getDims() const noexcept {
        if constexpr (Dim == 1)
            return IndexType{(unsigned char)superSize[0], numUnitCellSite};
        if constexpr (Dim == 2)
            return IndexType{(unsigned char)superSize[0], (unsigned char)superSize[1], numUnitCellSite};
        if constexpr (Dim == 3)
            return IndexType{(unsigned char)superSize[0], (unsigned char)superSize[1], (unsigned char)superSize[2], numUnitCellSite};
    }

    template<class Derived>
    bool LatticeHamilton<Derived>::checkNumSite() const noexcept {
        unsigned int numSite = getNumUnitCellSite();
        for (auto size : superSize)
            numSite *= size;
        return numSite == NumSite;
    }
}

namespace Physica {
    template<class Derived>
    class Traits<Core::LatticeHamilton<Derived>> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColumnAtCompile = Dynamic;
        constexpr static size_t MaxRowAtCompile = Dynamic;
        constexpr static size_t MaxColumnAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static size_t MaxSizeAtCompile = Dynamic;
    };
}
