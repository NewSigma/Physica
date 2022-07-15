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
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    template<class MeshType>
    class AbstractModel {
    public:
        using ScalarType = typename MeshType::ScalarType;
    protected:
        using ElementType = typename MeshType::ElementType;
        using VectorType = typename ElementType::VectorType;

        MeshType mesh;
        FEMSolver<ScalarType> solver;
    private:
        Utils::Array<size_t> map_var_node;
        Utils::Array<size_t> map_node_var;
    public:
        AbstractModel(MeshType mesh_);
        /* Operators */
        [[nodiscard]] ScalarType operator()(VectorType p) const { return mesh(p); }
        /* Getters */
        [[nodiscard]] const MeshType& getMesh() const noexcept { return mesh; }
    protected:
        /* Operations */
        void updateMesh();
        /* Getters */
        [[nodiscard]] size_t varToNode(size_t x) const { return map_var_node[x]; }
        [[nodiscard]] size_t nodeToVar(size_t node) const;
    private:
        /* Operations */
        void makeMaps();
    };

    template<class MeshType>
    AbstractModel<MeshType>::AbstractModel(MeshType mesh_) : mesh(std::move(mesh_)), solver() {
        const size_t n = mesh.getNumFreeNodes();
        solver.resize(n);
        map_var_node.resize(n);
        map_node_var.resize(mesh.getNumNodes());
        makeMaps();
    }

    template<class MeshType>
    size_t AbstractModel<MeshType>::nodeToVar(size_t node) const {
        assert(mesh.getNodeTypes()[node] == NodeType::Free);
        return map_node_var[node];
    }

    template<class MeshType>
    void AbstractModel<MeshType>::makeMaps() {
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
    void AbstractModel<MeshType>::updateMesh() {
        auto& coeffs = mesh.getCoeffs();
        for (size_t i = 0; i < map_var_node.getLength(); ++i) {
            coeffs[map_var_node[i]] = solver.x[i];
        }
    }
}
