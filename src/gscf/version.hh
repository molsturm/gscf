//
// Copyright (C) 2016 by the linalgwrap authors
//
// This file is part of linalgwrap.
//
// linalgwrap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// linalgwrap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "gscf/version_defs.hh"  // will be created by cmake
#include <string>

// The file gscf/version_defs.hh will be created by cmake from
// the interal project version and will contain definitions of the
// macros VERSION_MINOR VERSION_MAJOR VERSION_PATCH

namespace gscf {

struct version {
    static int constexpr major{gscf_VERSION_MAJOR};
    static int constexpr minor{gscf_VERSION_MINOR};
    static int constexpr patch{gscf_VERSION_PATCH};

    // Return the version as a string
    static std::string version_string();
};

}  // namespace gscf
