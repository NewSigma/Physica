/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "AbstractModel.h"

namespace Physica {
    template<class MeshType, class Functor>
    class PoissonModel : public AbstractModel<MeshType> {
    public:
        using Base = AbstractModel<MeshType>;
        using typename Base::ScalarType;
        using typename Base::ElementType;
    private:
        using Base::mesh;
        using typename Base::VectorType;

        Functor func;
    public:
        PoissonModel(MeshType mesh, Functor func_);
        /* Operations */
        template<class Integrator>
        void solve();
    };
    /**
     * \param func
     * ScalarType Functor(VectorType)
     */
    template<class MeshType, class Functor>
    PoissonModel<MeshType, Functor>::PoissonModel(MeshType mesh, Functor func_)
            : Base(std::move(mesh)), func(std::move(func_)) {}

    template<class MeshType, class Functor>
    template<class Integrator>
    void PoissonModel<MeshType, Functor>::solve() {
        Base::A.zeros();
        Base::b = ScalarType(0);
        const auto& nodeTypes = mesh.getNodeTypes();
        const auto& coeffs = mesh.getCoeffs();
        for (const auto& elem : mesh.getElements()) {
            const auto& nodes = elem.getGlobalNodes();
            for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                const size_t node = nodes[i];
                const bool isFree = nodeTypes[node] == NodeType::Free;
                if (isFree) {
                    const size_t row = Base::nodeToVar(node);

                    for (size_t j = 0; j < ElementType::getNumNodes(); ++j) {
                        const size_t baseNode = nodes[j];
                        const auto func = [=, &elem](VectorType p) {
                            const auto inv_jacobi = elem.inv_jacobi(p);
                            const VectorType g1 = inv_jacobi.transpose() * elem.grad(i, p);
                            const VectorType g2 = inv_jacobi.transpose() * elem.grad(j, p);
                            return abs(elem.jacobi(p).det()) * (g1 * g2);
                        };
                        const ScalarType integral = ElementType::template gauss_integral<-1>(func);

                        switch (nodeTypes[baseNode]) {
                        case NodeType::Free: {
                            const size_t col = Base::nodeToVar(baseNode);
                            const ScalarType value = Base::A.calc(row, col) + integral;
                            if (!value.isZero())
                                Base::A.insert(value, row, col);
                            break;
                        }
                        case NodeType::Dirichlet:
                            Base::b[row] -= coeffs[baseNode] * integral;
                            break;
                        }
                    }

                    Base::b[row] -= Integrator::run([&, i](VectorType p) {
                        return abs(elem.jacobi(p).det()) * elem.baseFunc(i, p) * func(elem.toGlobalPos(p));
                    });
                }
            }
        }
        Base::solve();
        Base::solverToMesh();
    }
}
