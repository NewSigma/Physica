/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include <memory>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunction.h"

namespace Physica::Core {
    /*!
     * @class AbstractIntegrate contains public members between integrate classes.
     */
    template<ScalarType type, bool errorTrack>
    class AbstractIntegrate {
    protected:
        std::shared_ptr<TreeFunction<type, errorTrack>> p_func;

        explicit AbstractIntegrate(std::shared_ptr<TreeFunction<type, errorTrack>>&& p);
    public:
        /* Getters */
        [[nodiscard]] const TreeFunction<type, errorTrack>& getFunction() const noexcept { return *p_func; }
    };

    template<size_t dim, ScalarType type = MultiPrecision, bool errorTrack = true>
    class Integrate : public AbstractIntegrate<type, errorTrack> {
        //We consider dimension which is larger than 3 is unsuitable to allocate domain data on stack.
        static_assert(dim <= 3, "Dimension larger than 3 must be set Dynamic.");
    };

    template<ScalarType type, bool errorTrack>
    class Integrate<Dynamic, type, errorTrack> : public AbstractIntegrate<type, errorTrack> {
        typedef AbstractIntegrate<type, errorTrack> Base;
        CStyleArray<Scalar<type, errorTrack>, Dynamic> vector;
    };

    template<ScalarType type, bool errorTrack>
    class Integrate<1, type, errorTrack> : public AbstractIntegrate<type, errorTrack> {
        typedef AbstractIntegrate<type, errorTrack> Base;
        Scalar<type, errorTrack> from, to;
    public:
        Integrate(std::shared_ptr<TreeFunction<type, errorTrack>> p
                  , Scalar<type, errorTrack> from, Scalar<type, errorTrack> to);
        /* Getters */
        const Scalar<type, errorTrack>& getFrom() const { return from; }
        const Scalar<type, errorTrack>& getTo() const { return to; }
    };

    template<ScalarType type, bool errorTrack>
    class Integrate<2, type, errorTrack> : public AbstractIntegrate<type, errorTrack> {
        typedef AbstractIntegrate<type, errorTrack> Base;
        Scalar<type, errorTrack> from1, to1, from2, to2;
    };

    template<ScalarType type, bool errorTrack>
    class Integrate<3, type, errorTrack> : public AbstractIntegrate<type, errorTrack> {
        typedef AbstractIntegrate<type, errorTrack> Base;
        Scalar<type, errorTrack> from1, to1, from2, to2, from3, to3;
    };
}

#include "IntegrateImpl.h"

#endif
