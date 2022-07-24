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
        using BufferType = typename Base::SolverType::VectorType;
        using Base::mesh;
        using Base::solver;

        Functor force;
        Functor diffuse;
        Functor d_diffuse;
        ScalarType stepSize;
        ScalarType mass;
        BufferType buffer1;
    public:
        FokkerPlanckModel(MeshType mesh,
                          Functor func_,
                          Functor diffuse_,
                          Functor d_diffuse_,
                          ScalarType stepSize_,
                          ScalarType mass_);
        /* Operations */
        void setInitialCond(Functor initial);
        void step();
    private:
        void spaceStep();
        void meshToBuffer(BufferType& buffer);
        void bufferToMesh(const BufferType& buffer);
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
            , mass(mass_) {
        buffer1.resize(Base::getDegreeOfFreedom());
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::setInitialCond(Functor initial) {
        solver.A.clear();
        solver.b = ScalarType::Zero();
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
                                const ScalarType value = solver.A.calc(row, col) + integral;
                                solver.A.insert(value, row, col);
                                break;
                            }
                            case NodeType::Dirichlet: {
                                solver.b[row] -= coeffs[baseNode] * integral;
                                break;
                            }
                        }
                    }

                    solver.b[row] += ElementType::gauss_integral([&, i](VectorType p) -> ScalarType {
                                         const VectorType globalPos = elem.toGlobalPos(p);
                                         const VectorType phase{tan(globalPos[0]), tan(globalPos[1])};
                                         if (std::isfinite(double(phase[0])) && std::isfinite(double(phase[1])))
                                            return abs(elem.jacobi(p).determinate()) * elem.baseFunc(i, p) * initial(phase);
                                         else
                                            return 0;
                                     });
                }
            }
        }
        solver.solve();
        Base::solverToMesh();
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::step() {
        meshToBuffer(buffer1);
        spaceStep();
        buffer1 += solver.b * stepSize;
        bufferToMesh(buffer1);
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::spaceStep() {
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
                        const ScalarType integral = ElementType::gauss_integral(
                                [=, &elem](VectorType p) {
                                    const VectorType globalPos = elem.toGlobalPos(p);
                                    const auto inv_jacobi = elem.inv_jacobi(p);
                                    const VectorType global_grad_i = inv_jacobi.transpose() * elem.grad(i, p);
                                    const VectorType global_grad_j = inv_jacobi.transpose() * elem.grad(j, p);

                                    const ScalarType cos_eta = cos(globalPos[1]);
                                    const bool flag = cos_eta < std::numeric_limits<ScalarType>::epsilon(); //Avoid divide by zero
                                    if (flag)
                                        return ScalarType::Zero();

                                    const ScalarType cos_xi = cos(globalPos[0]);
                                    const ScalarType sin_eta = sin(globalPos[1]);
                                    const ScalarType tan_xi = sin(globalPos[0]) / cos_xi;
                                    const ScalarType tan_eta = sin_eta / cos_eta;
                                    const VectorType phase_pos{tan_xi, tan_eta};

                                    const ScalarType u_i = elem.baseFunc(i, p);
                                    const ScalarType result1 = -tan_eta / mass * square(cos_xi) * u_i * global_grad_j[0];

                                    const ScalarType diffuseD = diffuse(phase_pos);
                                    const ScalarType d_diffuseD = d_diffuse(phase_pos);
                                    const ScalarType squared_cos_eta = square(cos_eta);
                                    const ScalarType temp = -cos_eta * sin_eta * u_i * 2 + squared_cos_eta * global_grad_i[1];
                                    const ScalarType result2 = (force(phase_pos) + (diffuseD - 1) * d_diffuseD) * (temp * elem.baseFunc(j, p));

                                    const ScalarType result3 = -(diffuseD * squared_cos_eta) * (temp * global_grad_j[1]);
                                    return (result1 + result2 + result3) * abs(elem.jacobi(p).determinate());
                                });
                        solver.b[row] += integral * mesh.getCoeffs()[baseNode];
                    }
                }
            }
        }
        solver.solve();
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::meshToBuffer(BufferType& buffer) {
        auto& coeffs = mesh.getCoeffs();
        for (size_t i = 0; i < Base::getDegreeOfFreedom(); ++i) {
            const size_t index = Base::varToNode(i);
            buffer[i] = coeffs[index];
        }
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::bufferToMesh(const BufferType& buffer) {
        auto& coeffs = mesh.getCoeffs();
        for (size_t i = 0; i < Base::getDegreeOfFreedom(); ++i) {
            const size_t index = Base::varToNode(i);
            coeffs[index] = buffer[i];
        }
    }
}
