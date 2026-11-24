/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cassert>
#include <utility>
#include <H5Ipublic.h>
#include "Physica/Macro.h"

namespace Physica {
    /**
     * \class H5ID: A type-safe RAII wrapper for hid_t of HDF5
     */
    class PHYSICA_API H5ID {
        hid_t id = H5I_INVALID_HID;
    public:
        enum class IdentifierType : int8_t {
            Invalid,
            File,
            Group,
            Datatype,
            Dataspace,
            Dataset,
            Attribute
        };
    public:
        H5ID() = default;
        explicit H5ID(hid_t id_) noexcept;
        H5ID(const H5ID& other) noexcept;
        H5ID(H5ID&& other) noexcept;
        ~H5ID() noexcept;
        /* Operators */
        H5ID& operator=(H5ID other) noexcept;
        /* Operations */
        template<class T>
        [[nodiscard]] bool isa() const noexcept;
        template<class T>
        [[nodiscard]] T cast(this H5ID&& self) noexcept;

        void checkWrite() const;
        void swap(H5ID& other) noexcept;
        /* Getters */
        [[nodiscard]] auto getHID() const noexcept { return id; }
        [[nodiscard]] bool isValid() const noexcept;

        [[nodiscard]] bool isReadOnly() const noexcept;
    protected:
        void incRef() const noexcept;
    private:
        [[nodiscard]] IdentifierType itype() const noexcept;
    };

    template<class T>
    bool H5ID::isa() const noexcept {
        return itype() == T::itype();
    }

    template<class T>
    T H5ID::cast(this H5ID&& self) noexcept {
        assert(self.isa<T>());
        return T(std::move(self));
    }
}
