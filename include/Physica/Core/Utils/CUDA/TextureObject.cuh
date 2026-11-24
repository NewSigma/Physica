/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Macro.h"

namespace Physica {
    template<class T>
    class TextureObject {
        using This = TextureObject;

        cudaTextureObject_t texObj;
    public:
        TextureObject() : texObj(0) {}
        TextureObject(size_t major, size_t minor);
        TextureObject(const This&) = delete;
        TextureObject(This&& other) noexcept;
        ~TextureObject();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ operator cudaTextureObject_t() const { return texObj; }
        /* Operations */
        void fromMatrix(const Matrix auto& m);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] cudaResourceDesc getResourceDesc() const noexcept;
        [[nodiscard]] cudaArray_t getArray() const noexcept { return getResourceDesc().res.array.array; }
    };

    template<class T>
    TextureObject<T>::TextureObject(size_t major, size_t minor) {
        cudaResourceDesc resDesc;
        /* Make resDesc */ {
            memset(&resDesc, 0, sizeof(resDesc));
            resDesc.resType = cudaResourceTypeArray;

            cudaChannelFormatDesc channelDesc;
            if constexpr (std::is_same<T, float16>::value)
                channelDesc = cudaCreateChannelDesc<int16_t>(); // cudaCreateChannelDescHalf() seems do not work using CUDA 12.6
            else if constexpr (std::is_same<T, float32>::value)
                channelDesc = cudaCreateChannelDesc<float>();
            else if constexpr (std::is_same<T, float64>::value)
                channelDesc = cudaCreateChannelDesc<double>();
            else
                channelDesc = cudaCreateChannelDesc<T>();
            check(cudaMallocArray(&resDesc.res.array.array, &channelDesc, minor, major));
        }

        cudaTextureDesc texDesc;
        /* Make texDesc */ {
            memset(&texDesc, 0, sizeof(texDesc));
            texDesc.addressMode[0] = cudaAddressModeClamp;
            texDesc.addressMode[1] = cudaAddressModeClamp;
            texDesc.filterMode = cudaFilterModePoint;
            texDesc.readMode = cudaReadModeElementType;
            texDesc.normalizedCoords = 0;
        }
        check(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));
        assert(texObj != 0 && "[Error]: We have assumed null is 0, but it is not the case");
    }

    template<class T>
    TextureObject<T>::TextureObject(This&& other) noexcept : texObj(other.texObj) {
        other.texObj = 0;
    }

    template<class T>
    TextureObject<T>::~TextureObject() {
        if (texObj != 0) {
            cudaFreeArray(getArray());
            cudaDestroyTextureObject(texObj);
        }
    }

    template<class T>
    void TextureObject<T>::fromMatrix(const Matrix auto& m) {
        const size_t major = m.getMaxMajor();
        const size_t minor = m.getMaxMinor();
        const size_t spitch = minor * sizeof(T);
        check(cudaMemcpy2DToArray(getArray(), 0, 0, m.data(), spitch, spitch, major, cudaMemcpyHostToDevice));
    }

    template<class T>
    void TextureObject<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(texObj, obj.texObj);
    }

    template<class T>
    cudaResourceDesc TextureObject<T>::getResourceDesc() const noexcept {
        cudaResourceDesc result;
        cudaGetTextureObjectResourceDesc(&result, texObj);
        return result;
    }
}
