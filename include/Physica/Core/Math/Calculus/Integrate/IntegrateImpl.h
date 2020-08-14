/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_INTEGRATEIMPL_H
#define PHYSICA_INTEGRATEIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    AbstractIntegrate<type, errorTrack>::AbstractIntegrate(std::shared_ptr<TreeFunction<type, errorTrack>>&& p)
            : p_func(std::move(p)) {}
    ////////////////////////////////////Dynamic////////////////////////////////////
    ////////////////////////////////////1D////////////////////////////////////
    template<ScalarType type, bool errorTrack>
    Integrate<1, type, errorTrack>::Integrate(std::shared_ptr<TreeFunction<type, errorTrack>> p
            , Scalar<type, errorTrack> from, Scalar<type, errorTrack> to)
                : AbstractIntegrate<type, errorTrack>(std::move(p)), from(std::move(from)), to(std::move(to)) {}
    ////////////////////////////////////2D////////////////////////////////////
    ////////////////////////////////////3D////////////////////////////////////
}

#endif
