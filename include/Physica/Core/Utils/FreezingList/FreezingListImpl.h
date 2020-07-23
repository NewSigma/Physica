/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FREEZINGLISTIMPL_H
#define PHYSICA_FREEZINGLISTIMPL_H
/*!
 * This File is the implementation of FreezingList,
 * Do not include this file, include FreezingList.h instead.
 */
namespace Physica::Core {
    /*!
     * Allocate the space, but no data is stored.
     *
     * Warning: Do not attempt to access the data.
     */
    template<class T>
    FreezingList<T>::FreezingList(size_t length) : arr(reinterpret_cast<Node*>(malloc(length * sizeof(Node))))
            , length(length), inited(0) {
        start = arr;
    }

    template<class T>
    FreezingList<T>::~FreezingList() {
        free(arr);
    }
    /*!
     * Allocate data in sequence.
     */
    template<class T>
    void FreezingList<T>::allocate(const T& t) {
        const auto copy = inited;
        if(copy < length) {
            ++inited;
            Node* p = arr + copy;
            if(QTypeInfo<T>::isComplex)
                new (p) T(t);
            else
                p->t = t;
            arr[copy].last = (copy == 0) ? nullptr : p - 1;
            arr[copy].next = (inited == length) ? nullptr : p + 1;
        }
    }

    template<class T>
    void FreezingList<T>::allocate(T&& t) {
        const auto copy = inited;
        if(copy < length) {
            ++inited;
            Node* p = arr + copy;
            if(QTypeInfo<T>::isComplex)
                new (p) T(std::move(t));
            else
                p->t = t;
            arr[copy].last = (copy == 0) ? nullptr : p - 1;
            arr[copy].next = (inited == length) ? nullptr : p + 1;
        }
    }
    /*!
     * Freeze the node if the node belongs to this list.
     */
    template<class T>
    void FreezingList<T>::freeze(Node* node) {
        if(node >= arr && node - arr < length) {
            Node* n_last = node->last;
            Node* n_next = node->next;
            start = start == node ? n_next : start;
            if(n_last)
                n_last->next = n_next == nullptr ? nullptr : n_next;
            if(n_next)
                n_next->last = n_last == nullptr ? nullptr : n_last;
        }
    }
    /*!
     * Restore the node if the node belongs to this list.
     */
    template<class T>
    void FreezingList<T>::restore(Node* node) {
        if(node >= arr && node - arr < length) {
            //List is empty.
            if(start == nullptr) {
                start = node;
                node->last = node->next = nullptr;
                return;
            }
            //node is before front.
            if(node < start) {
                node->last = nullptr;
                node->next = start;
                start->last = node;
                start = node;
                return;
            }
            //node is between front and the end of the list.
            Node* p = start;
            Node* p_next = p->next;
            while(p_next != nullptr) {
                if(p < node && node < p_next) {
                    node->last = p;
                    node->next = p_next;
                    p_next->last = node;
                    p->next = node;
                    return;
                }
                p = p_next;
                p_next = p->next;
            }
            //node is after the end of the list.
            p->next = node;
            node->last = p;
            node->next = nullptr;
        }
    }
}

#endif