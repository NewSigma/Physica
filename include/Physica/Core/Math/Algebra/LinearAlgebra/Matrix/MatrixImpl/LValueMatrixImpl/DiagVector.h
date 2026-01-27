/*
 * Copyright 2025-2026 Weibo He.
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

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class DiagVectorL : public LValueVector<DiagVectorL<M>> {
        using This = DiagVectorL<M>;
        using Base = LValueVector<This>;
    private:
        LazyDestroy<M> mat;
    public:
        explicit DiagVectorL(M&& mat) : mat(std::forward<M>(mat)) {}
        DiagVectorL(const This&) = default;
        DiagVectorL(This&&) = default;
        ~DiagVectorL() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t index) noexcept { return self.mat.data_ptr(index, index); }
    };

    template<Matrix M>
    auto&& DiagVectorL<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<DiagVectorL<M>> : public Traits<DiagVectorR<M>> {};
}
