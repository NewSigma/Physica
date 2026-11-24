/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Calculus/PDE/FEM/Mesh.h"

namespace Physica {
    template<class ElementType>
    class VTKFile {
        using MeshType = Mesh<ElementType>;

        const MeshType& mesh;
        std::string title;
    public:
        VTKFile(const MeshType& mesh_, std::string title_);
        /* Operators */
        template<class T>
        friend std::ostream& operator<<(std::ostream& os, const VTKFile<T>& vtk);
    };

    template<class ElementType>
    VTKFile<ElementType>::VTKFile(const MeshType& mesh_, std::string title_) : mesh(mesh_), title(std::move(title_)) {}

    template<class ElementType>
    std::ostream& operator<<(std::ostream& os, const VTKFile<ElementType>& vtk) {
        const auto& mesh = vtk.mesh;
        const size_t numGlobalNode = mesh.getNumNodes();
        const size_t numElems = mesh.getNumElems();
        /* Title */ {
            os << "# vtk DataFile Version 3.0\n"
               << vtk.title << '\n'
               << "ASCII\n";
        }
        /* Geometry */ {
            os << "DATASET POLYDATA\n"
               << "POINTS " << numGlobalNode << " float\n";
            const auto nodes = mesh.getNodes();
            for (const auto& node : nodes) {
                constexpr const char* suffix = ElementType::Dim == 2 ? " 0\n" : "\n";
                os << node.format().setPrefix("").setSeparator(" ").setSuffix(suffix);
            }

            constexpr size_t nodePerElem = ElementType::getNumNodes();
            const size_t arraySize = (nodePerElem + 1) * numElems;
            os << "POLYGONS " << numElems << ' ' << arraySize << '\n';
            for (const auto& element : mesh.getElements()) {
                os << nodePerElem;
                for (size_t node : element.getGlobalNodes())
                    os << ' ' << node;
                os << '\n';
            }
        }
        /* Data */ {
            const auto& coeffs = mesh.getCoeffs();
            os << "POINT_DATA " << numGlobalNode << '\n';
            os << "SCALARS fieldName double 1\n"
               << "LOOKUP_TABLE default\n"
               << coeffs.format().setPrefix("").setSeparator(" ").setSuffix("\n") << '\n';
        }
        return os;
    }
}
