#ifndef PHYSICA_NET_H
#define PHYSICA_NET_H

#include <vector>
#include <unordered_map>
#include "Layer.h"
#include "Vector.h"

class Net {
private:
    class Connections {
        struct interface {
            int node;
            int interface;
        };
    private:
        Net* belongs;
        typedef std::unordered_multimap<int, interface> connection;
        int size;
        std::unordered_map<int, connection>** connections;
    public:
        Connections(Net* belongs, int size);
        ~Connections();

        bool connect(int from_layer, int from_node, int to_layer, int to_node, int interface);
        void disconnect(int from_layer, int from_node, int to_layer, int to_node, int interface);
    };

    Vector* variables;
    Layer** layers;
    Connections* connections;
    int length;
public:
    Net(Vector* variables, std::vector<int> structure);
    ~Net();

    Layer* getLayer(int index);
};

#endif