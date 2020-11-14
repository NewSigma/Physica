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
    namespace Utils {
        /*!
         * CUUniquePtr is the unique_ptr in CUDA,
         * which is allocated on host and handles pointers on GPU.
         *
         * T must be a trivial type because we can not call destructor functions on device.
         */
        template<class T>
        class CUUniquePtr {
            T* pointer;
        public:
            CUUniquePtr(T* t) : pointer(t) {}
            CUUniquePtr(const CUUniquePtr<T>&) = delete;
            CUUniquePtr(CUUniquePtr<T>&& ptr) : pointer(ptr.pointer) { ptr.pointer = 0; }
            ~CUUniquePtr() { cudaFree(pointer); }
            /* Operators */
            CUUniquePtr& operator=(const CUUniquePtr<T>) = delete;
            CUUniquePtr& operator=(CUUniquePtr<T>&& ptr);
            T& operator*() { return *pointer; }
            T* operator->() { return pointer; }
            /* Operations */
            T* get() { return pointer; }
            T* release();
            void reset(T* p);
            void swap(CUUniquePtr<T>& ptr);
        };

        template<class T>
        CUUniquePtr<T>& CUUniquePtr<T>::operator=(CUUniquePtr<T>&& ptr) {
            ~CUUniquePtr();
            pointer = ptr.pointer;
            ptr.pointer = 0;
        }

        template<class T>
        T* CUUniquePtr<T>::release() {
            auto temp = pointer;
            pointer = 0;
            return temp;
        }

        template<class T>
        void CUUniquePtr<T>::reset(T* p) {
            auto temp = pointer;
            pointer = p;
            p = temp;
            if(p != pointer)
                ~CUUniquePtr();
        }

        template<class T>
        void CUUniquePtr<T>::swap(CUUniquePtr<T>& ptr) {
            auto temp = pointer;
            pointer = ptr.pointer;
            ptr.pointer = temp;
        }
    }
}

#endif
