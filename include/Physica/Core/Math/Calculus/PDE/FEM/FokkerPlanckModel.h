/*
 * Copyright 2022-2024 Weibo He.
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
        using BufferType = VectorND<ScalarType>;
        using Base::mesh;

        Functor force;
        Functor diffuse;
        Functor d_diffuse;
        ScalarType stepSize;
        ScalarType rep_mass;
        SparseMatrix<ScalarType> stiff;
        BufferType buffer1;
        BufferType buffer2;
    public:
        FokkerPlanckModel(MeshType mesh,
                          Functor func_,
                          Functor diffuse_,
                          Functor d_diffuse_,
                          ScalarType stepSize_,
                          ScalarType mass);
        /* Operations */
        void setInitialCond(Functor initial);
        void step();
    private:
        void makeStiff();
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
                                                            ScalarType mass)
            : Base(std::move(mesh))
            , force(std::move(force_))
            , diffuse(std::move(diffuse_))
            , d_diffuse(std::move(d_diffuse_))
            , stepSize(stepSize_)
            , rep_mass(reciprocal(mass)) {
        const size_t n = Base::getDegreeOfFreedom();
        stiff = SparseMatrix<ScalarType>(n, n);
        buffer1.resize(Base::getDegreeOfFreedom());
        buffer2.resize(Base::getDegreeOfFreedom());
        makeStiff();
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::setInitialCond(Functor initial) {
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
                        const ScalarType integral = ElementType::gauss_integral(
                                [=, &elem](VectorType p) {
                                    return abs(elem.jacobi(p).determinate()) * elem.baseFunc(i, p) * elem.baseFunc(j, p);
                                });

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

                    const auto func = [&, i](VectorType p) -> ScalarType {
                        const VectorType globalPos = elem.toGlobalPos(p);
                        return abs(elem.jacobi(p).determinate()) * elem.baseFunc(i, p) * initial(globalPos);
                    };
                    Base::b[row] += ElementType::template gauss_integral<decltype(func), 0>(func);
                }
            }
        }
        Base::solve();
        Base::solverToMesh();
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::step() {
        meshToBuffer(buffer1);
        buffer2 = buffer1;
        spaceStep();
        buffer1 += Base::b * (stepSize * 0.5);
        spaceStep();
        buffer2 += Base::b * stepSize;
        bufferToMesh(buffer2);
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::makeStiff() {
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
                        const auto func = [=, this, &elem](VectorType p) {
                            const VectorType globalPos = elem.toGlobalPos(p);
                            const auto inv_jacobi = elem.inv_jacobi(p);
                            const VectorType global_grad_i = inv_jacobi.transpose() * elem.grad(i, p);
                            const VectorType global_grad_j = inv_jacobi.transpose() * elem.grad(j, p);

                            const ScalarType u_i = elem.baseFunc(i, p);
                            const ScalarType result1 = -globalPos[1] * rep_mass * u_i * global_grad_j[0];

                            const ScalarType diffuseD = diffuse(globalPos);
                            const ScalarType d_diffuseD = d_diffuse(globalPos);
                            const ScalarType result2 = (force(globalPos) + (diffuseD - 1) * d_diffuseD) * (elem.baseFunc(j, p) * global_grad_i[1]);

                            const ScalarType result3 = -diffuseD * (global_grad_i[1] * global_grad_j[1]);
                            return (result1 + result2 + result3) * abs(elem.jacobi(p).determinate());
                        };
                        const ScalarType integral = ElementType::template gauss_integral<decltype(func), 0>(func);
                        const size_t col = Base::nodeToVar(baseNode);
                        const ScalarType value = stiff.calc(row, col) + integral;
                        stiff.insert(value, row, col);
                    }
                }
            }
        }
    }

    template<class MeshType, class Functor>
    void FokkerPlanckModel<MeshType, Functor>::spaceStep() {
        Base::b = stiff * buffer1;
        Base::solve();
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
