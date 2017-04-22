//
// Copyright (C) 2017 by the gscf authors
//
// This file is part of gscf.
//
// gscf is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// gscf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with gscf. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "gscf/version_defs.hh"  // will be created by cmake
#include <string>

// The file gscf/version_defs.hh will be created by cmake from
// the interal project version and will contain definitions of the
// macros VERSION_MINOR VERSION_MAJOR VERSION_PATCH

namespace gscf {

struct version {
  static int constexpr major{detail::__version_var_major};
  static int constexpr minor{detail::__version_var_minor};
  static int constexpr patch{detail::__version_var_patch};

  // Return the version as a string
  static std::string version_string();
};

}  // namespace gscf
