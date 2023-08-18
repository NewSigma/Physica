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

namespace Physica::Core {
    template<class SpaceType, size_t Dim> class SubDataSpace;

    namespace Internal {
        template<class T, size_t S>
        class Traits<SubDataSpace<T, S>> {
        public:
            constexpr static size_t Dim = S;
        };
    }

    template<class SpaceType, size_t Dim>
    class SubDataSpace : public DataSpaceBase<SubDataSpace<SpaceType, Dim>> {
        static_assert(Dim != Utils::Dynamic, "[Error]: Not implemented");
        using This = SubDataSpace<SpaceType, Dim>;
        using Base = DataSpaceBase<This>;
    public:
        using SizeArray = typename Base::SizeArray;
    private:
        SpaceType& space;
        size_t fromDim;
        size_t toDim;
    public:
        SubDataSpace(DataSpaceBase<SpaceType>& space_, size_t fromDim_, size_t toDim_);
        SubDataSpace(const SubDataSpace&) = default;
        SubDataSpace(SubDataSpace&&) noexcept = default;
        ~SubDataSpace() = default;
        /* Operations */
        void selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start);
        /* Getters */
        [[nodiscard]] const H5::DataSpace& asH5Type() const noexcept { return space.asH5Type(); }
        [[nodiscard]] H5::DataSpace& asH5Type() noexcept { return space.asH5Type(); }
        [[nodiscard]] inline size_t getDim() const noexcept;
        [[nodiscard]] size_t getSize(size_t dim) const noexcept { return space.getSize(fromDim + dim); }
    };

    template<class SpaceType, size_t Dim>
    SubDataSpace<SpaceType, Dim>::SubDataSpace(DataSpaceBase<SpaceType>& space_, size_t fromDim_, size_t toDim_)
            : space(space_.getDerived())
            , fromDim(fromDim_)
            , toDim(toDim_) {
        assert((Dim == Dynamic || (Dim == toDim - fromDim)) && "[Error]: Inconsistent between type and its param");
    }

    template<class SpaceType, size_t Dim>
    void SubDataSpace<SpaceType, Dim>::selectHyperslab(H5S_seloper_t op, const SizeArray& count, const SizeArray& start) {
        using ExtendSizeArray = typename SpaceType::SizeArray;
        assert((op == H5S_seloper_t::H5S_SELECT_SET) && "[Error]: Not implemented");
        ExtendSizeArray count1 = space.getSelectedCount();
        ExtendSizeArray start1 = space.getSelectedStart();
        for (size_t i = 0; i < getDim(); ++i) {
            count1[i + fromDim] = count[i];
            start1[i + fromDim] = start[i];
        }
        space.selectHyperslab(op, count1, start1);
    }

    template<class SpaceType, size_t Dim>
    inline size_t SubDataSpace<SpaceType, Dim>::getDim() const noexcept {
        if constexpr (Dim == Utils::Dynamic)
            return toDim - fromDim;
        else
            return Dim;
    }
}
