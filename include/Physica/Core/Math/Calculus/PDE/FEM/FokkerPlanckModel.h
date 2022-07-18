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
    class FokkerPlanckModel : public AbstractModel<MeshType> {
    public:
        using Base = AbstractModel<MeshType>;
        using typename Base::ScalarType;
        using typename Base::ElementType;
    private:
        using typename Base::VectorType;
        using Base::mesh;
        using Base::solver;

        Functor force;
        Functor diffuse;
        Functor d_diffuse;
        ScalarType stepSize;
        ScalarType mass;
    public:
        FokkerPlanckModel(MeshType mesh,
                          Functor func_,
                          Functor diffuse_,
                          Functor d_diffuse_,
                          ScalarType stepSize_,
                          ScalarType mass_);
        /* Operations */
        template<class Integrator>
        void setInitialCond(Functor initial);
        template<class Integrator>
        void step();
    private:
        void solverToMesh();
        void updateMesh();
    };

    template<class MeshType, class Functor>
    FokkerPlanckModel<MeshType, Functor>::FokkerPlanckModel(MeshType mesh,
                                                            Functor force_,
                                                            Functor diffuse_,
                                                            Functor d_diffuse_,
                                                            ScalarType stepSize_,
                                                            ScalarType mass_)
            : Base(std::move(mesh))
            , force(std::move(force_))
            , diffuse(std::move(diffuse_))
            , d_diffuse(std::move(d_diffuse_))
            , stepSize(stepSize_)
            , mass(mass_) {}

    template<class MeshType, class Functor>
    template<class Integrator>
    void FokkerPlanckModel<MeshType, Functor>::setInitialCond(Functor initial) {
        solver.clear();
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
                        const ScalarType integral = ElementType::gauss_integral(
                                [=, &elem](VectorType p) {
                                    return abs(elem.jacobi(p).determinate()) * elem.baseFunc(i, p) * elem.baseFunc(j, p);
                                });

                        switch (nodeTypes[baseNode]) {
                            case NodeType::Free: {
                                const size_t col = Base::nodeToVar(baseNode);
                                solver.A(row, col) += integral;
                                break;
                            }
                            case NodeType::Dirichlet: {
                                solver.b[row] -= coeffs[baseNode] * integral;
                                break;
                            }
                        }
                    }

                    solver.b[row] += Integrator::run([&, i](VectorType p) {
                                         return abs(elem.jacobi(p).determinate()) * elem.baseFunc(i, p) * initial(elem.toGlobalPos(p));
                                     });
                }
            }
        }
        solver.solve();
        for (auto& xi : solver.x)
            if (xi.isNegative())
                xi = ScalarType::Zero();
        Base::solverToMesh();
    }

    template<class MeshType, class Functor>
    template<class Integrator>
    void FokkerPlanckModel<MeshType, Functor>::step() {
        solver.b = ScalarType::Zero();
        const auto& nodeTypes = mesh.getNodeTypes();
        for (const auto& elem : mesh.getElements()) {
            const auto& nodes = elem.getGlobalNodes();
            for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                const size_t node = nodes[i];
                const bool isFree = nodeTypes[node] == NodeType::Free;
                if (isFree) {
                    const size_t row = Base::nodeToVar(node);
                    for (size_t j = 0; j < ElementType::getNumNodes(); ++j) {
                        const size_t baseNode = nodes[j];
                        const ScalarType integral = Integrator::run(
                                [=, &elem](VectorType p) {
                                    const VectorType globalPos = elem.toGlobalPos(p);
                                    const auto inv_jacobi = elem.inv_jacobi(p);
                                    const VectorType global_grad_i = inv_jacobi.transpose() * elem.grad(i, p);
                                    const VectorType global_grad_j = inv_jacobi.transpose() * elem.grad(j, p);

                                    const ScalarType term1 = ScalarType::One() - square(globalPos[0]);
                                    const ScalarType term2 = ScalarType::One() - square(globalPos[1]);
                                    const ScalarType term3 = term1 * sqrt(term1);

                                    const bool flag = term2 < std::numeric_limits<ScalarType>::epsilon(); //Avoid divide by zero
                                    const ScalarType u_i = elem.baseFunc(i, p);
                                    const ScalarType result1 = flag ? ScalarType::Zero()
                                                                    : (-globalPos[1] / term2 * term3 * (u_i * global_grad_j[0]));

                                    const ScalarType diffuseD = diffuse(globalPos);
                                    const ScalarType d_diffuseD = d_diffuse(globalPos);
                                    const ScalarType result2 = (force(globalPos) + (diffuseD - 1) * d_diffuseD) * sqrt(term2) * ((term2 * global_grad_i[1] - ScalarType(3) * globalPos[1] * u_i) * elem.baseFunc(j, p));

                                    const ScalarType squaredTerm2 = square(term2);
                                    const ScalarType result3 = ScalarType(3) * diffuseD * (globalPos[1] * squaredTerm2) * (u_i * global_grad_j[1]);
                                    const ScalarType result4 = -diffuseD * squaredTerm2 * term2 * (global_grad_i[1] * global_grad_j[1]);
                                    return (result1 + result2 + result3 + result4) * abs(elem.jacobi(p).determinate());
                                });
                        solver.b[row] += integral * mesh.getCoeffs()[baseNode];
                    }
                }
            }
        }
        solver.solve();
        updateMesh();
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::updateMesh() {
        auto& coeffs = mesh.getCoeffs();
        for (size_t i = 0; i < Base::getDegreeOfFreedom(); ++i) {
            const size_t index = Base::varToNode(i);
            coeffs[index] += solver.x[i] * stepSize;
        }
    }
}
