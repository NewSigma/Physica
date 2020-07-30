/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTCSTYLEARRAYIMPL_H
#define PHYSICA_ABSTRACTCSTYLEARRAYIMPL_H

namespace Physica::Core {
    template<class T>
    typename AbstractCStyleArray<T>::Iterator& AbstractCStyleArray<T>::Iterator::operator=(const Iterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class T>
    typename AbstractCStyleArray<T>::Iterator& AbstractCStyleArray<T>::Iterator::operator=(Iterator&& ite) noexcept {
        p = ite.p;
        return *this;
    }

    template<class T>
    typename AbstractCStyleArray<T>::Iterator& AbstractCStyleArray<T>::Iterator::operator++() {
        ++p;
        return *this;
    }

    template<class T>
    const typename AbstractCStyleArray<T>::Iterator AbstractCStyleArray<T>::Iterator::operator++(int) { //NOLINT Return type must be const.
        return Iterator(p++);
    }

    template<class T>
    AbstractCStyleArray<T>::AbstractCStyleArray(std::initializer_list<T> list)
            : arr(reinterpret_cast<T*>(malloc(list.size() * sizeof(T)))), length(list.size()) {
        size_t i = 0;
        const auto end = list.end();
        for(auto ite = list.begin(); ite != end; ++ite, ++i)
            allocate(T(*ite), i);
    }

    template<class T>
    AbstractCStyleArray<T>::AbstractCStyleArray(AbstractCStyleArray<T>&& array) noexcept
            : arr(array.arr), length(array.length) {
        array.arr = nullptr;
        array.length = 0;
    }

    template<class T>
    AbstractCStyleArray<T>::~AbstractCStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T>
    AbstractCStyleArray<T>& AbstractCStyleArray<T>::operator=(AbstractCStyleArray<T>&& array) noexcept {
        this->~AbstractCStyleArray();
        arr = array.arr;
        length = array.length;
        array.arr = nullptr;
        array.length = 0;
        return *this;
    }

    template<class T>
    bool AbstractCStyleArray<T>::operator==(const AbstractCStyleArray& array) const {
        if(length != array.length)
            return false;
        for(size_t i = 0; i < length; ++i)
            if(operator[](i) != array[i])
                return false;
        return true;
    }
    /*!
     * Low level api. Designed for performance.
     * Simply allocate a \T at position \index.
     * You must ensure position \index is usable or a memory leak will occur.
     */
    template<class T>
    inline void AbstractCStyleArray<T>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void AbstractCStyleArray<T>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(std::move(t));
        else
            *(arr + index) = t;
    }
    /*!
     * Get the last element in the array and remove it from the array.
     */
    template<class T>
    T AbstractCStyleArray<T>::cutLast() {
        Q_ASSERT(length > 0);
        --length;
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }

    template<class T>
    void AbstractCStyleArray<T>::swap(AbstractCStyleArray<T> &array) {
        std::swap(arr, array.arr);
        std::swap(length, array.length);
    }
}

#endif
