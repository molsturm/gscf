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
#include "ScfBaseKeys.hh"

namespace gscf {

/** Struct which contains the keys used for setting the
 *  TruncatedOptDampScf parameters */
struct TruncatedOptDampScfKeys : public ScfBaseKeys {
  /** The number of previous SCF steps to consider for the truncated optimal
   *  damping SCF algorithm.
   *  Type: size_type */
  static const std::string n_prev_steps;

  /** The smallest prefactor which may sit in front of a Fock matrix
   *  for it to be considered at all (Type: scalar_type) */
  static const std::string min_fock_prefactor;
};

}  // namespace gscf
