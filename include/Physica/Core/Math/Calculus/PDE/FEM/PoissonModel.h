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

#include "FEMSolver.h"
#include "Mesh.h"

namespace Physica::Core {
    template<class MeshType>
    class PoissonModel : public FEMSolver<typename MeshType::ScalarType> {
    public:
        using ScalarType = typename MeshType::ScalarType;
        using ElementType = typename MeshType::ElementType;
    private:
        using Base = FEMSolver<ScalarType>;
        using VectorType = typename ElementType::VectorType;

        MeshType mesh;
        Utils::Array<size_t> map_var_node;
        Utils::Array<size_t> map_node_var;
    public:
        PoissonModel(MeshType mesh_);
        /* Operators */
        [[nodiscard]] ScalarType operator()(VectorType p) const { return mesh(p); }
        /* Operations */
        template<class Functor, class Integrator>
        void solve(Functor func);
        /* Getters */
        [[nodiscard]] const MeshType& getMesh() const noexcept { return mesh; }
    private:
        /* Operations */
        void makeMaps();
        void postSolve();
        /* Getters */
        [[nodiscard]] size_t varToNode(size_t x) const { return map_var_node[x]; }
        [[nodiscard]] size_t nodeToVar(size_t node) const;
    };

    template<class MeshType>
    PoissonModel<MeshType>::PoissonModel(MeshType mesh_) : Base(), mesh(std::move(mesh_)) {
        const size_t n = mesh.getNumFreeNodes();
        Base::resize(n);
        map_var_node.resize(n);
        map_node_var.resize(mesh.getNumNodes());
        makeMaps();
    }
    /**
     * \param func
     * ScalarType Functor(VectorType)
     */
    template<class MeshType>
    template<class Functor, class Integrator>
    void PoissonModel<MeshType>::solve(Functor func) {
        Base::clear();

        const auto& nodeTypes = mesh.getNodeTypes();
        const auto& coeffs = mesh.getCoeffs();
        for (const auto& elem : mesh.getElements()) {
            const auto& nodes = elem.getGlobalNodes();
            for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                const size_t node = nodes[i];
                const bool isFree = nodeTypes[node] == NodeType::Free;
                if (isFree) {
                    const size_t row = nodeToVar(node);

                    for (size_t j = 0; j < ElementType::getNumNodes(); ++j) {
                        const size_t baseNode = nodes[j];
                        const ScalarType integral = ElementType::integral(
                                [=, &elem](VectorType p) {
                                    const auto inv_jacobi = elem.inv_jacobi(p);
                                    const VectorType g1 = inv_jacobi.transpose() * elem.grad(i, p);
                                    const VectorType g2 = inv_jacobi.transpose() * elem.grad(j, p);
                                    return abs(elem.jacobi(p).determinate()) * (g1 * g2);
                                });

                        switch (nodeTypes[baseNode]) {
                            case NodeType::Free: {
                                const size_t col = nodeToVar(baseNode);
                                Base::A(row, col) += integral;
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
        postSolve();
    }

    template<class MeshType>
    void PoissonModel<MeshType>::makeMaps() {
        const auto& nodeTypes = mesh.getNodeTypes();
        size_t next_x = 0;
        for (size_t i = 0; i < nodeTypes.getLength(); ++i) {
            const bool isFreeNode = nodeTypes[i] == NodeType::Free;
            if (isFreeNode) {
                map_var_node[next_x] = i;
                map_node_var[i] = next_x;
                ++next_x;
            }
        }
    }

    template<class MeshType>
    void PoissonModel<MeshType>::postSolve() {
        auto& coeffs = mesh.getCoeffs();
        for (size_t i = 0; i < map_var_node.getLength(); ++i) {
            coeffs[map_var_node[i]] = Base::x[i];
        }
    }

    template<class MeshType>
    size_t PoissonModel<MeshType>::nodeToVar(size_t node) const {
        assert(mesh.getNodeTypes()[node] == NodeType::Free);
        return map_node_var[node];
    }
}
