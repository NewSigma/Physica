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
    template<size_t Dim>
    class H5DataSet : public H5::DataSet {
        using Base = H5::DataSet;
        constexpr static auto Dynamic = Physica::Utils::Dynamic;
    public:
        H5DataSet() = default;
        H5DataSet(const H5::DataSet& obj);
        H5DataSet(const H5DataSet&) = default;
        H5DataSet(H5DataSet&&) noexcept = delete;
        virtual ~H5DataSet() = default;
        /* Operators */
        inline H5DataSet& operator=(const H5DataSet& obj);
        H5DataSet& operator=(H5DataSet&&) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] H5DataSpace<Dim> getDataSpace() const noexcept { return Base::getSpace(); }
        [[nodiscard]] size_t getDim() const noexcept;
        [[nodiscard]] size_t getSize(size_t dim) const noexcept;
    private:
        using Base::getSpace;
    };

    template<size_t Dim>
    H5DataSet<Dim>::H5DataSet(const H5::DataSet& obj) : Base(obj) {}

    template<size_t Dim>
    inline H5DataSet<Dim>& H5DataSet<Dim>::operator=(const H5DataSet& obj) {
        Base::operator=(obj);
        return *this;
    }

    template<size_t Dim>
    size_t H5DataSet<Dim>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return getDataSpace().getDim();
        else
            return Dim;
    }

    template<size_t Dim>
    size_t H5DataSet<Dim>::getSize(size_t dim) const noexcept {
        return getDataSpace().getSize(dim);
    }
}
