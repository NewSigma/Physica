/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_CU_UNIQUE_PTR_H
#define PHYSICA_CU_UNIQUE_PTR_H

namespace Physica {
    namespace Core {
        /*!
         * CUUniquePtr is the unique_ptr in CUDA.
         */
        template<class T>
        class CUUniquePtr {
            T* pointer;
        public:
            __device__ CUUniquePtr(T* t) : pointer(t) {}
            __device__ CUUniquePtr(const CUUniquePtr<T>&) = delete;
            __device__ CUUniquePtr(CUUniquePtr<T>&& ptr) : pointer(ptr.pointer) { ptr.pointer = 0; }
            __device__ ~CUUniquePtr() { delete pointer; }
            /* Operators */
            __device__ CUUniquePtr& operator=(const CUUniquePtr<T>) = delete;
            __device__ CUUniquePtr& operator=(CUUniquePtr<T>&& ptr);
            __device__ T& operator*() { return *pointer; }
            __device__ T* operator->() { return pointer; }
            /* Operations */
            __device__ T* get() { return pointer; }
            __device__ T* release();
            __device__ void reset(T* p);
            __device__ void swap(CUUniquePtr<T>& ptr);
        };

        template<class T>
        __device__ CUUniquePtr<T>& CUUniquePtr<T>::operator=(CUUniquePtr<T>&& ptr) {
            ~CUUniquePtr();
            pointer = ptr.pointer;
            ptr.pointer = 0;
        }

        template<class T>
        __device__ T* CUUniquePtr<T>::release() {
            auto temp = pointer;
            pointer = 0;
            return temp;
        }

        template<class T>
        __device__ void CUUniquePtr<T>::reset(T* p) {
            auto temp = pointer;
            pointer = p;
            p = temp;
            if(p != pointer)
                ~CUUniquePtr();
        }

        template<class T>
        __device__ void CUUniquePtr<T>::swap(CUUniquePtr<T>& ptr) {
            auto temp = pointer;
            pointer = ptr.pointer;
            ptr.pointer = temp;
        }
    }
}

#endif
