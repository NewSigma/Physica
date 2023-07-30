/*
 * Copyright 2022 WeiBo He.
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

#include "AbstractModel.h"

namespace Physica::Core {
    template<class MeshType, class Functor>
    class PoissonModel : public AbstractModel<MeshType> {
    public:
        using Base = AbstractModel<MeshType>;
        using typename Base::ScalarType;
        using typename Base::ElementType;
    private:
        using typename Base::VectorType;
        using Base::mesh;

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
    PoissonModel<MeshType, Functor>::PoissonModel(MeshType mesh, Functor func_) : Base(std::move(mesh)), func(std::move(func_)) {}

    template<class MeshType, class Functor>
    template<class Integrator>
    void PoissonModel<MeshType, Functor>::solve() {
        Base::A.clear();
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
                            return abs(elem.jacobi(p).determinate()) * (g1 * g2);
                        };
                        const ScalarType integral = ElementType::template gauss_integral<decltype(func), -1>(func);

                        switch (nodeTypes[baseNode]) {
                            case NodeType::Free: {
                                const size_t col = Base::nodeToVar(baseNode);
                                const ScalarType value = Base::A.calc(row, col) + integral;
                                Base::A.insert(value, row, col);
                                break;
                            }
                            case NodeType::Dirichlet: {
                                Base::b[row] -= coeffs[baseNode] * integral;
                                break;
                            }
                        }
                    }

                    Base::b[row] -= Integrator::run([&, i](VectorType p) {
                                        return abs(elem.jacobi(p).determinate()) * elem.baseFunc(i, p) * func(elem.toGlobalPos(p));
                                    });
                }
            }
        }
        Base::solve();
        Base::solverToMesh();
    }
}
