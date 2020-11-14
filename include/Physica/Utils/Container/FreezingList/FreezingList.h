/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_FREEZINGLIST_H
#define PHYSICA_FREEZINGLIST_H

#include <cstddef>

namespace Physica::Utils {
    /*!
     * Comparision:
     * This class is similar to std::vector because it is stored consistently and
     * similar to std::list because it has list nodes.
     * The only difference is that FreezingList can freeze one of its nodes or restore a frozen node,
     * If we freeze the n.th node, the list behaves as if the node does not exist and the (n - 1).th node
     * will directly link to the (n + 1).th node.
     *
     * Same to CStyleArray, only classes that have copy or move constructor can be stored in this container.
     * This container never call the default constructor of T.
     *
     * Used in:
     * TotalCircuit
     */
    template<class T>
    class FreezingList {
    public:
        struct Node {
            T t;
            Node* last;
            Node* next;
        };

        class Iterator {
            Node* node;
        public:
            ~Iterator() = default;

            Iterator(const Iterator& ite) : node(ite.node) {}
            Iterator(Iterator&& ite) noexcept : node(ite.node) {}

            /* Operators */
            Iterator& operator=(const Iterator& ite);
            Iterator& operator=(Iterator&& ite) noexcept;
            bool operator==(const Iterator& ite) const noexcept { return node == ite.node; }
            bool operator!=(const Iterator& ite) const noexcept { return node != ite.node; }
            Iterator& operator++();
            const Iterator operator++(int);
            T& operator*() { return node->t; }
        private:
            explicit Iterator(Node* node) : node(node) {}
            Node* getNode() const { return node; }

            friend class FreezingList;
        };
    private:
        //Array of nodes, allocated by malloc.
        Node* arr;
        //The beginning of list, different from arr because arr may be frozen.
        Node* start;
        //Length of arr.
        size_t length;
        //Length of initialize arr.
        size_t inited;
    public:
        explicit FreezingList(size_t length);
        FreezingList(std::initializer_list<T> list);
        FreezingList(const FreezingList& list) = delete;
        FreezingList(FreezingList&& list) noexcept = delete;
        ~FreezingList();
        /* Operators */
        FreezingList& operator=(const FreezingList& list) = delete;
        FreezingList& operator=(FreezingList&& list) noexcept = delete ;
        /* Operations */
        void allocate(const T& t);
        void allocate(T&& t);
        void freeze(const Iterator& ite);
        void restore(const Iterator& ite);
        /* Getters */
        //Data must be initialized or the access will be refused.
        [[nodiscard]] size_t size() const { return length; }
        [[nodiscard]] bool empty() const { return start == nullptr; }
        /* Iterator */
        [[nodiscard]] Iterator begin() const { Iterator(this->start); }
        [[nodiscard]] Iterator end() const { return Iterator(nullptr); }
    };
}

#include "FreezingListImpl.h"

#endif