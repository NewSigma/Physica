/*
 * Copyright 2019 WeiBo He.
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
#ifndef PHYSICA_LAYER_H
#define PHYSICA_LAYER_H

#include <vector>

namespace Physica::AI {
    class DNN;

    class Layer {
        std::vector<Node*> nodes;
        DNN* parentNet;
        //index of current layer in the net
        int id;
    public:
        explicit Layer(int id, int nodeCount, DNN* parent);
        ~Layer();
        Layer(const Layer&) = delete;
        Layer& operator=(const Layer&) = delete;
        /* Operators */
        Node& operator[](int i) { return *nodes[i]; }
        /* Getters */
        [[nodiscard]] DNN* getNet() const { return parentNet; }
        [[nodiscard]] int getSize() const { return nodes.size(); }
        [[nodiscard]] int getId() const { return id; }
        /* Operations */
        void update();
        void handleLoss();
    };
}

#endif
